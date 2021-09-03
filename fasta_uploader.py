""" fasta_uploader.py

See latest at: https://github.com/Public-Health-Bioinformatics/FastaUploader/

Given a fasta file and a sample metadata file with a column that matches to
fasta file record identifiers, break both into respective sets of smaller
batches of records which are submitted to an API for processing.

Processing is three step: 

1) construct batches of files. Since two files are read and parsed in one go,
processing of them is reliable after that point, so no further error reporting
required after parsing.
 - Importantly, if rerunning, this step will be skipped unless -f --force 
 parameter is run. Currently input files are still required in this case.

2) IF API option is included, submit each batch to API, wait for it to finish
or error out (capture error report) and proceed to next batch. Successfully
uploaded batches are renamed to include their job ids by the API.
- Some types of error trigger sudden death, i.e. sys.exit() because they 
would apply to any subsequent API batch calls.  For example missing tabular
data column names will trigger an exit. Once resolved, rerun with -f to force
regeneration of output files.

3) Scan through all uploaded fasta batches and report back via the API any new
information about errors they may contain.

Authors: Damion Dooley, Nithu Sara John
Centre for Infectious Disease Epidemiology and One Health
August 24, 2021

Requires Biopython and Requests modules
- "pip install biopython" or https://biopython.org/wiki/Download
- "pip install requests"

Usage:
python fasta_uploader.py -c "20210713_AB_final set 1.csv" -f "consensus_renamed_final.fasta" -k "fasta header name" -d -a VirusSeq_Portal -u eyJhbGciOiJSUz....ykapL1A

FUTURE: Add feature to remerge all split fasta files and tsv files to enable
them to be error corrected in batch files directly, rather than correcting
problems in original merged file and rerunning whole upload process?
"""

from Bio import SeqIO
import numpy as np
import optparse
import pandas as pd
import sys
import requests
import os, glob
from datetime import datetime


# init_parser()
# returns command line parameters
#
# @return tuple options, args: Dict options command line parameters by name
#
def init_parser():

   parser = optparse.OptionParser();

   parser.add_option('-f', '--fasta', dest="fasta_file",
      help="Provide a fasta file name.");
   parser.add_option('-c', '--csv', dest="csv_file", 
      help="Provide a COMMA delimited sample contextual data file name.");
   parser.add_option('-t', '--tsv', dest="tsv_file", 
      help="Provide a TAB delimited sample contextual data file name.");
   parser.add_option('-b', '--batch', dest="batch",
      help="Provide number of fasta records to include in each batch.", default=1000);
   parser.add_option('-o', '--output', dest="output_file",
      help="Provide an output file name/path.", default='output');
   parser.add_option('-k', '--key', dest="key_field",
      help="Provide the metadata field name to match to fasta record identifier.");
   parser.add_option('-n', '--number', dest="batch_number", default = False,
      help="Just process given batch number to API instead of all batches.");
   # API related parameters

   parser.add_option('-a', '--api', dest="api",
      help="provide the target API to send data too.  A batch submission job will be initiated for it.");
   parser.add_option('-u', '--user', dest="api_token",
      help="An API user token is required for API access.");
   parser.add_option('-d', '--dev', dest="development", action='store_true',
      help="Test against a development server rather than live one.  Provide an API endpoint URL");
   parser.add_option('-s', '--short', dest="short", default = False, action='store_false',
      help="Report up to given # of fasta record related errors for each batch submission.  Useful for taking care of repeated errors first based on first instance.");
   parser.add_option('-r', '--reset', dest="reset", action='store_true', 
      help="Regenerate all batch files and begin API resubmission process even if batch files already exist under given output file pattern.");

   return parser.parse_args();


# log()
# Presents textual messages to standard i/o, but also writes them to a log file.
# 
#
def log(log_handler, text): 
   print (text);
   log_handler.write(text + '\n');
   return text;


# get_metadata()
# Reads a .tsv or .csv file containing sample contextual data and sorts it
#
# @param Object log_handler for saving progress and error text
# @param Dict options command line parameters by name
# @return list sorted by options.key_field
#
def get_metadata(log_handler, options):

   if not options.fasta_file:
      sys.exit(log(log_handler, "A sample sequencing fasta file is required."));

   if options.csv_file:
      metadata = pd.read_csv(options.csv_file, encoding = 'unicode_escape');
      file_suffix = 'csv';

   elif options.tsv_file:
      metadata = pd.read_table(options.tsv_file, delimiter='\t', encoding = 'unicode_escape');
      file_suffix = 'tsv';

   else:
      sys.exit(log(log_handler, "A sample contextual .tsv or .csv file is required."));

   # Check if given fasta record identifier is a sample metadata file column
   if not options.key_field in metadata.columns:
      log(log_handler, 'The key field column you provided [' + options.key_field + '] was not found in the contextual data file\'s list of columns:');
      log(log_handler, metadata.columns);
      sys.exit(1);

   metadata.sort_values(by = options.key_field);

   return metadata;


# get_fasta_data()
# Creates sorted list (array) of fasta records.
#
# @param Object log_handler for saving progress and error text
# @param Dict options command line parameters by name
#
def get_fasta_data(log_handler, options):

   with open(options.fasta_file, "r") as fasta_handle:
      fasta_data = SeqIO.parse(fasta_handle, "fasta");
      # Sort Fasta file so we organize upload, and can sync with metadata
      fasta_data = [f for f in sorted(fasta_data, key=lambda x : x.id)];

   return fasta_data;


# batch_fasta()
# Creates [options.output_file].X records containing fasta and sample 
# contextual data in batches of at most [options.batch] records each.
#
# @param Object log_handler for saving progress and error text
# @param Object fasta_data SeqIO object containing all fasta records
# @param List metadata list containing all fasta contextual data records
# @param Dict options command line parameters by name
#
def batch_fasta(log_handler, fasta_data, metadata, options):
   
   # Splits into batches of options.batch (default 1000) or less records:
   splits = len(fasta_data)/int(options.batch);
   if splits < 1:
      splits = 1;

   for count, sequences in enumerate(np.array_split(fasta_data, splits)):
      
      # Determine metadata rows pertinent to all sequences. They should be in same order
      write_metadata(log_handler, sequences, metadata, count, options);

      with open(options.output_file + '.'+ str(count) + '.queued.fasta', 'w') as output_handle:
         SeqIO.write(sequences, output_handle, "fasta");


# write_metadata
#
# @param Object log_handler for saving progress and error text

# @param Object fasta_data SeqIO object containing all fasta records
# @param List metadata list containing all fasta contextual data records
# @param String count contains batch file 
# @param Dict options command line parameters by name
# @param String id last parameter of tsv or csv file name, 'queued' by default 
#
def write_metadata(log_handler, fasta_data, metadata, count, options, id='queued'):
   id_index = [record.id for record in fasta_data];

   metabatch = metadata.loc[metadata[options.key_field].isin(id_index)];

   #log(log_handler, 'Files for ' + options.output_file + '.'+ str(count))
   # https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.to_csv.html
   if options.csv_file:
      metabatch.to_csv(options.output_file + '.'+ str(count) + '.'+id+'.csv', index=False); 
   else:
      # Or to tabular?  Issue of quoted strings?
      metabatch.to_csv(options.output_file + '.'+ str(count) + '.'+id+'.tsv', sep="\t", index=False);


# api_batch_job()
# Submits two files, fasta and sample contextual data, to given API endpoint 
# at a time for all files in current working directory.
#
# @param Object log_handler for saving progress and error text
# @param Object fasta_data SeqIO object containing all fasta records
# @param List metadata list containing all fasta contextual data records
# @param Dict options command line parameters by name
#
def api_batch_job(log_handler, fasta_data, metadata, options):

   if not options.api_token:
      sys.exit(log(log_handler, "An API user token is required for use with the [" + options.api + "] API."));

   # Retrieve batch files that need uploading, indicated by ".queued." in filename.
   batches = glob.glob("./" + options.output_file + '.*.queued.fasta');

   if len(batches) > 0:
      log(log_handler, 'Performing batch upload ...');
      
      # Currently only virusseq () API is an option.
      if options.api == 'VirusSeq_Portal':

         api_virusseq_job(log_handler, options, batches);

   else:
      log(log_handler, 'Skipping API batch file queue because it is empty.');


################################## VirusSeq API ##########################
#
#   API Information: https://github.com/cancogen-virus-seq/docs/wiki/How-to-Submit-Data-(API)
#
#   Log in here to get API Key, good for 24 hours:
#   https://portal.dev.cancogen.cancercollaboratory.org/login
#
#   NOTE: Study_id field if not set to account project like "UHTC-ON" will
#   trigger validation error "UNAUTHORIZED_FOR_STUDY_UPLOAD".
#

def api_virusseq_job(log_handler, options, batches):

   if options.development:  # TEST API ENDPOINT
      API_URL = "https://muse.dev.cancogen.cancercollaboratory.org/";
   else:  # LIVE API ENDPOINT
      API_URL = "https://muse.virusseq-dataportal.ca/";  

   # TEST file parse failure: create an empty or junky .fasta and accompanying .tsv file
   # batches = ['test.0.queued.fasta'];

   for filename in batches:

      count = filename.split('.')[-3];

      if (not options.batch_number or count == options.batch_number):

         filename_tsv = filename.replace('.queued.fasta', '.queued.tsv');
         with open(filename, 'rb') as fasta_handle:
            with open(filename_tsv, 'rb') as metadata_handle:
               upload_files = [
                  ('files', fasta_handle), 
                  ('files', metadata_handle)
               ];
               log(log_handler, 'Uploading batch: ' + filename);
               try:
                  request = requests.post(API_URL + 'submissions', files = upload_files, headers = {'Authorization': 'Bearer ' + options.api_token});
               except Exception as err:
                  sys.exit(log(log_handler, "API Server problem (check API URL?): " + repr(err) ));

         if request.status_code == 200:
            result = request.json();
            if ('submissionId' in result):
               submission_id = result['submissionId'];
               log(log_handler, 'Batch was submitted! submissionId: ' + submission_id);

               os.rename(filename, filename.replace('.queued.','.'+submission_id+'.'));
               os.rename(filename_tsv, filename_tsv.replace('.queued.','.'+submission_id+'.'));

               continue;
            else:
               log(log_handler, result);
               sys.exit(log(log_handler, "Resolve reported error, then rerun command!"));

         # "Unauthorized client error status response" code occurs when key is not valid.
         # This halts processing of all remaining batches.
         if request.status_code == 401:
            log(log_handler, "Unauthorized client error status 401 response");
            sys.exit(log(log_handler, "Check to make sure your API key is current."));

         if request.status_code == 404:
            sys.exit(log(log_handler, "API service endpoint not recognized. Check API URL:" + API_URL))
         
         request_error = request.json();
         status = request_error['status'];
         message = request_error['message'];

         log(log_handler, request_error);
         errorInfo = request_error['errorInfo'];

         # "Bad Request" response status indicates something wrong with the input files.
         # We have json at this point.
         if request.status_code == 400: # or status == "BAD_REQUEST"
            log(log_handler, "Bad Request status 400 response.");

            if (status == "BAD_REQUEST"):

               if message == 'Headers are incorrect!':
                  log(log_handler, message);
                  log(log_handler, "Unknown Headers:", errorInfo['unknownHeaders']);
                  log(log_handler, "Missing Headers:", errorInfo['missingHeaders']);
                  sys.exit(log(log_handler, "Check to make sure the .tsv file headers are current."));

               """
               Example error:
               "errorInfo":{"invalidFields":[{"fieldName":"specimen collector sample ID","value":"","reason":"NOT_ALLOWED_TO_BE_EMPTY","index":1},{"fieldName":"fasta header name","value":"","reason":"NOT_ALLOWED_TO_BE_EMPTY","index":1},{"fieldName":"study_id","value":" 23434","reason":"UNAUTHORIZED_FOR_STUDY_UPLOAD","index":1}]}}
               """
               if message == 'Found records with invalid fields':
                  for record in errorInfo['invalidFields']:
                     log(log_handler, "row " + str(record['index']), '"' + record['fieldName'] + '"', 
                        record['reason'],
                        "value:",   record['value']
                     );
                  continue;

            # not sure where this should be positioned.
            # {"status":"FORBIDDEN","message":"Denied","errorInfo":{}}
            if (status == "FORBIDDEN"):
               sys.exit(log(log_handler, "Your account associated with the API key has not been authorized, so this service is not available to you yet."));

         # Internal Server Error (code generated etc.)
         if request.status_code == 500:
            log(log_handler, status);
            log(log_handler, message);
            if message == "Flux#last() didn't observe any onNext signal":
               log(log_handler, "Does .tsv file have no data rows?");
            continue;

         log(log_handler, 'Error: Unable to complete batch because of status code ' + str(request.status_code) + '\n' + request.text);
         continue;


def api_batch_status(log_handler, fasta_data, metadata, options):
   
   # Get list of batch files to get status for
   batches = glob.glob('./' + options.output_file + '.*.*.fasta');
   #batches = batches.sort(key=os.path.getmtime);
   batches = sorted(batches, key=lambda n: int(n.split('.',3)[2]) ); # sort on count of file
   for filename in batches:
      
      # Any errors get converted into entries in the following files
      error_keys = [];

      file_name_parts = filename.split('.');
      submission_id = filename.split('.')[-2];
      count = filename.split('.')[-3];

      if (not options.batch_number or count == options.batch_number):

         if not submission_id == 'queued':
            log(log_handler, '\nSTATUS for: ' + filename);
            if options.short:
               error_max = options.short;
            else:
               error_max = options.batch;

            if (options.api == 'VirusSeq_Portal'):

               query = '?page=0&size=' + str(error_max) + '&sortDirection=ASC&sortField=submitterSampleId&submissionId=' + submission_id;

               if options.development:  # TEST API ENDPOINT
                  API_URL = "https://muse.dev.cancogen.cancercollaboratory.org/";
               else:  # LIVE API ENDPOINT
                  API_URL = "https://muse.virusseq-dataportal.ca/";  

               feedback = requests.get(API_URL + 'uploads' + query, headers = {'Authorization': 'Bearer ' + options.api_token});

               if feedback.status_code == 200:
                  response = feedback.json();

                  item_report = '';

                  for submission in response['data']:

                     # We wait for results of queued job item
                     if (submission['status'] == 'QUEUED'):
                        item_report += submission['submitterSampleId'] + " Queued" + '\n';

                     # Any error item must be resubmitted
                     if (submission['status'] == 'ERROR'):
                        error_key = submission['submitterSampleId'];
                        error_keys.append(error_key);

                        for ptr, item in enumerate(submission['error'].split('#')[1:]):
                           # For brevity, just show field name, not section
                           binding = item.split(':',1);
                           item_label = binding[0].split('/')[-1];
                           item_error = binding[1];
                           item_report += '\n' + error_key + '\t' + item_label + item_error;

                  log(log_handler, item_report);

               else:
                  log(log_handler, 'Status unavailable');

               log(log_handler, '\n');

         # Write out appropriate fasta, metadata and file content for each error.
         # This causes only errors to be resubmitted - AFTER EDITING - on next run.

         if len(error_keys):
            # A FILTER ITERATOR
            sequences = list(filter(lambda x: x.id in error_keys, fasta_data));

            # np.array(fasta_data)[in error_keys]:
            with open(options.output_file + '.'+ count + '.queued.fasta', 'w') as output_handle:

               SeqIO.write(sequences, output_handle, "fasta");
               # just getting header, no content ????
               # sequences is EMPTY?
               write_metadata(log_handler, sequences, metadata, count, options);
            
            # TRY: Update job files to only include successes
            sequences = list(filter(lambda x: not x.id in error_keys, fasta_data));
            with open(filename, 'w') as output_handle:

               SeqIO.write(sequences, output_handle, "fasta");
               write_metadata(log_handler, sequences, metadata, count, options, submission_id);

# ISSUE: NEW TSV GENERATION NOT BASED ON SUBMITTED FILE BUT RATHER ORIGINAL .tsv file.
# FIX by loading from [EDITED] submitted .tsv

################################## The Program ##########################

options, args = init_parser(); 

# Program operates within log file handler that closes on exit.
simple_date = datetime.now().isoformat().replace(':','').split('.')[0];
with open(options.output_file + '_' + simple_date + '.log', 'w') as log_handler:

   # Caution: Big file data loaded into memory
   fasta_data = get_fasta_data(log_handler, options);
   metadata = get_metadata(log_handler, options);

   log(log_handler, "Log date: " + datetime.now().isoformat() );

   # If force, clear out all existing batch files and start fresh
   if options.reset:
      # Add syntax check / security on options.output_file references?
      for filename in glob.glob("./" + options.output_file + '.*'):
         os.remove(filename);

   # STEP 1: GENERATE BATCH FILES
   # The existence of options.output_file + .*.*.fasta files present
   # indicates an ongoing job th
   batches = glob.glob("./" + options.output_file + '.0.*.fasta');
   if len(batches) > 0:
      log(log_handler, 'Skipping batch file generation because batch files exist.');
   else:
      log(log_handler, 'Generating batch file(s) ...' );
      batch_fasta(log_handler, fasta_data, metadata, options);

   # STEP 2: SUBMIT TO API
   if options.api:
      api_batch_job(log_handler, fasta_data, metadata, options);

   # STEP 3: Report on progress of each batch job that has been submitted.
      api_batch_status(log_handler, fasta_data, metadata, options);

