""" fasta_uploader.py

See latest at: https://github.com/cidgoh/FastaUploader/

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
information about errors they may contain.  For each batch file produce a new
.queued.fasta and .queued.tsv file containing records that have errors in them
so that they can be corrected manually, and then resubmitted by running
program again.

Authors: Damion Dooley, Nithu Sara John
Centre for Infectious Disease Epidemiology and One Health
August 24, 2021

Requires Biopython and Requests modules
- "pip install biopython" or https://biopython.org/wiki/Download
- "pip install requests"

Usage:
python fasta_uploader.py -m "samples_xyz.csv" -f "samples_xyz.fasta" -k "fasta header name" -d -a VirusSeq_Portal -u eyJhbGciOiJSUz....ykapL1A

FUTURE: Add feature to remerge all split fasta files and tsv files to enable
them so after error correction in batch files they can be merged into one
corrected file?
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
   parser.add_option('-m', '--metadata', dest="metadata_file", 
      help="Provide a COMMA .csv or TAB .tsv delimited sample contextual data file name.");
   parser.add_option('-b', '--batch', dest="batch",
      help="Provide number of fasta records to include in each batch.", default=1000);
   parser.add_option('-o', '--output', dest="output_file",
      help="Provide an output file name/path.", default='output');
   parser.add_option('-k', '--key', dest="key_field",
      help="Provide the metadata field name to match to fasta record identifier.");
   parser.add_option('-n', '--number', dest="batch_number", default = False,
      help="Process only given batch number to API instead of all batches.");
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
def get_metadata(log_handler, metadata_file, options):

   if not options.fasta_file:
      log(log_handler, "A sample sequencing fasta file is required.");
      sys.exit();

   if metadata_file[-4:].lower() == '.csv':
      metadata = pd.read_csv(metadata_file, encoding = 'unicode_escape');

   elif metadata_file[-4:].lower() == '.tsv':
      metadata = pd.read_table(metadata_file, delimiter='\t', encoding = 'unicode_escape');

   else:
      log(log_handler, "A sample contextual .tsv or .csv file is required.")
      sys.exit(1);

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
def get_fasta_data(log_handler, fasta_file, options):

   with open(fasta_file, "r") as fasta_handle:
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
   
   # Resort fasta and metadata so they have same order as split files.
   with open(options.fasta_file, 'w') as source_handle:
      log(log_handler, 'Sorting and resaving source fasta and tabular files');
      SeqIO.write(fasta_data, source_handle, "fasta");

      if options.metadata_file[-4:] == '.csv':
         separator = ',';
      else:
         separator = "\t";
      
      metadata.to_csv(options.metadata_file, sep=separator, index=False);


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
# see https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.to_csv.html
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

   # FOR NOW, ONLY .TSV OUTPUT for VirusSeq.
   #suffix = options.metadata_file[-4:];
   suffix = '.tsv';

   file_name = options.output_file + '.'+ str(count) + '.' + id + suffix;
   log(log_handler, 'Saving ' + file_name);

   if suffix == '.csv':
      separator = ',';
   else:
      separator = "\t";
   
   # Queued file might get several related batch number job status fail items added
   if os.path.exists(file_name) and id == 'queued':
      metabatch.to_csv(file_name, sep=separator, index=False, mode = 'a', header = False);
   else:
      metabatch.to_csv(file_name, sep=separator, index=False);


# api_batch_job()
# Submits two files, fasta and sample contextual data, to given API endpoint 
# at a time for all files in current working directory.
#
# @param Object log_handler for saving progress and error text
# @param Dict options command line parameters by name
#
def api_batch_job(log_handler, options):

   if not options.api_token:
      log(log_handler, "An API user token is required for use with the [" + options.api + "] API.");
      sys.exit(1);

   # Retrieve batch files that need uploading, indicated by ".queued." in filename.
   batches = glob.glob("./" + options.output_file + '.*.queued.fasta');

   if len(batches) > 0:
      log(log_handler, 'Performing batch upload ...');
      
      # Currently only virusseq () API is an option.
      if options.api == 'VirusSeq_Portal':

         api_virusseq_job(log_handler, options, batches);

   else:
      log(log_handler, 'API batch file queue is empty.');


################################## VirusSeq API ##########################
#
#   API Information: https://github.com/cancogen-virus-seq/docs/wiki/How-to-Submit-Data-(API)
#
#   Log in here to get API Key, good for 24 hours:
#   https://portal.dev.cancogen.cancercollaboratory.org/login
#
#   NOTE: Study_id field if not set to account project like UHTC-ON or DRGN-CA will
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
                  log(log_handler, "API Server problem (check API URL?): " + repr(err) );
                  sys.exit(1);

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
               log(log_handler, "Resolve reported error, then rerun command!");
               sys.exit(1);

         # "Unauthorized client error status response" code occurs when key is not valid.
         # This halts processing of all remaining batches.
         if request.status_code == 401:
            log(log_handler, "Unauthorized client error status 401 response");
            log(log_handler, "Check to make sure your API key is current.");
            sys.exit(1);

         if request.status_code == 404:
            log(log_handler, "API service endpoint not recognized. Check API URL:" + API_URL)
            sys.exit(1);
         
         request_error = request.json();
         status = request_error['status'];
         message = request_error['message'];

         log(log_handler, status + ' (' + str(request.status_code) + ') ' + message);

         errorInfo = request_error['errorInfo'];

         # "Bad Request" response status indicates something wrong with the input files.
         # We have json at this point.
         if request.status_code == 400: # or status == "BAD_REQUEST"
            log(log_handler, "Bad Request status 400 response.");

            if (status == "BAD_REQUEST"):

               if message == 'Headers are incorrect!':
                  log(log_handler, message);
                  log(log_handler, "Unknown Headers: " + str(errorInfo['unknownHeaders']));
                  log(log_handler, "Missing Headers: " + str(errorInfo['missingHeaders']));
                  log(log_handler, "Check to make sure the .tsv file headers are current.");
                  sys.exit(1);

               """
               Example error:
               "errorInfo":{"invalidFields":[{"fieldName":"specimen collector sample ID","value":"","reason":"NOT_ALLOWED_TO_BE_EMPTY","index":1},{"fieldName":"fasta header name","value":"","reason":"NOT_ALLOWED_TO_BE_EMPTY","index":1},{"fieldName":"study_id","value":" 23434","reason":"UNAUTHORIZED_FOR_STUDY_UPLOAD","index":1}]}}
               """
               if message == 'Found records with invalid fields':
                  for record in errorInfo['invalidFields']:
                     log(log_handler, "row " + str(record['index']) + ' "' + record['fieldName'] + '" ' + 
                        record['reason'] +
                        " value: " + record['value']);
                  continue;

            # not sure where this should be positioned.
            # {"status":"FORBIDDEN","message":"Denied","errorInfo":{}}
            if (status == "FORBIDDEN"):
               log(log_handler, "Your account associated with the API key has not been authorized, so this service is not available to you yet.")
               sys.exit(1);

         # Internal Server Error (code generated etc.)
         if request.status_code == 500:
            log(log_handler, status);
            log(log_handler, message);
            if message == "Flux#last() didn't observe any onNext signal":
               log(log_handler, "Does .tsv file have no data rows?");
            continue;

         log(log_handler, 'Error: Unable to complete batch because of status code ' + str(request.status_code) + '\n' + request.text);
         continue;


def api_batch_status(log_handler, options):
   
   # Get list of batch files to get status for
   batches = glob.glob('./' + options.output_file + '.*.*.fasta');
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

            #
            #if options.short:
            #   error_max = options.short;
            #else:
            #   error_max = options.batch;


            if (options.api == 'VirusSeq_Portal'):
               api_virusseq_status(log_handler, submission_id, error_keys, options)

         # Write out appropriate fasta and metadata for each error.
         # This causes only errors to be resubmitted - AFTER EDITING - on next run.

         # ISSUE: first row
         if len(error_keys) > 0:
            print (error_keys)

            fasta_data = get_fasta_data(log_handler, filename, options);
            filename_tsv = filename.replace('.fasta', '.tsv');
            metadata = get_metadata(log_handler, filename_tsv, options);

            # A filter iterator converted to a list
            sequences = list(filter(lambda x: x.id in error_keys, fasta_data));

            # Several batch id jobs can contribute to a single new .queued. file
            # if their own entry statuses switch from 'queued' to error
            output_fasta_file = options.output_file + '.'+ count + '.queued.fasta';
            if os.path.exists(output_fasta_file):
                append_write = 'a' # append if already exists
            else:
                append_write = 'w' # make a new file if not

            with open(output_fasta_file, append_write) as output_handle:

               SeqIO.write(sequences, output_handle, "fasta");
               write_metadata(log_handler, sequences, metadata, count, options);
            
            # Update job files to only include successes. Existing files are overwritten.
            sequences = list(filter(lambda x: not x.id in error_keys, fasta_data));
            with open(filename, 'w') as output_handle:

               SeqIO.write(sequences, output_handle, "fasta");
               write_metadata(log_handler, sequences, metadata, count, options, submission_id);


def api_virusseq_status(log_handler, submission_id, error_keys, options):

   # VirusSeq "size" parameter clips off last record status unless set to batch
   # size + 1.
   query = '?page=0&size=' + str(int(options.batch)+1) + '&sortDirection=ASC&sortField=submitterSampleId&submissionId=' + submission_id;

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
         else:
            if (submission['status'] == 'ERROR'):
               error_key = submission['submitterSampleId'];
               error_keys.append(error_key);

               for ptr, item in enumerate(submission['error'].split('#')[1:]):
                  # For brevity, just show field name, not section
                  binding = item.split(':',1);
                  item_label = binding[0].split('/')[-1];
                  item_error = binding[1];
                  item_report += '\n' + error_key + '\t' + item_label + item_error;

            else:
               item_report += submission['submitterSampleId'] + ' ' + submission['status'] + '\n';

      log(log_handler, item_report);

   else:
      log(log_handler, 'Status unavailable');

   log(log_handler, '\n');


################################## The Program ##########################

options, args = init_parser(); 

# Program operates within log file handler that closes on exit.
simple_date = datetime.now().isoformat().replace(':','').split('.')[0];
with open(options.output_file + '_' + simple_date + '.log', 'w') as log_handler:

   # Caution: Big file data loaded into memory
   fasta_data = get_fasta_data(log_handler, options.fasta_file, options);
   metadata = get_metadata(log_handler, options.metadata_file, options);

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

   if options.api:

      # STEP 2: SUBMIT the *.queued.* files TO API
      api_batch_job(log_handler, options);

      # STEP 3: Report on progress of each batch job that has been submitted.
      # Any previous error reports for a given batch job are prepared into a
      # new .queued.fasta and .queued.tsv file for editing and resubmission.
      api_batch_status(log_handler, options);

