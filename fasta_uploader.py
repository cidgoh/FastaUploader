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
python fasta_uploader.py -c "20210713_AB_final set 1.csv" -f "consensus_renamed_final.fasta" -k "fasta header name"

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

parser = optparse.OptionParser()

parser.add_option('-f', '--fasta', dest="fasta_file",
   help="provide a fasta file name");
parser.add_option('-c', '--csv', dest="csv_file", 
   help="provide a COMMA delimited sample contextual data file name");
parser.add_option('-t', '--tsv', dest="tsv_file", 
   help="provide a TAB delimited sample contextual data file name");
parser.add_option('-b', '--batch', dest="batch",
   help="provide number of fasta records to include in each batch", default=1000);
parser.add_option('-o', '--output', dest="output_file",
   help="provide an output file name/path", default='output');
parser.add_option('-k', '--key', dest="key_field",
   help="provide the metadata field name to match to fasta record identifier");

# API related parameters

parser.add_option('-a', '--api', dest="api", default='VirusSeq_Portal',
   help="provide the target API to send data too.  A batch submission job will be initiated for it.");
parser.add_option('-u', '--user', dest="api_token",
   help="an API user token is required for API access");
parser.add_option('-d', '--dev', dest="development",
   help="Test against a development server rather than live one.  Provide an API endpoint URL");
parser.add_option('-s', '--short', dest="short", default = False, 
   help="Report up to given # of fasta record related errors for each batch submission.  Useful for taking care of repeated errors first based on first instance.");
parser.add_option('-r', '--reset', dest="reset",
   help="regenerate all batch files and begin API resubmission process even if batch files already exist under given output file pattern.");

options, args = parser.parse_args();

with open(options.output_file + datetime.now().isoformat().replace(':','_') + '_log_.txt') as log_handle:

   metadata = get_metadata(log_handle, options);
   log(log_handle, "Fasta batch file process initiated on: " + datetime.now().isoformat() );

   # If force, clear out all existing batch files and redo
   if options.reset:
      # Add syntax check / security on options.output_file references?
      for filename in glob.glob("./" + options.output_file + '.*'):
         os.remove(filename);


   # STEP 1: GENERATE BATCH FILES
   batches = glob.glob("./" + options.output_file + '.*.*.fasta');

   # Doesn't regenerate fasta or tsv batches if any one matching pattern exists.
   if len(batches) > 0:
      print ('Skipping batch file generation because batch files exist.');

   else:


      print ('Generating batch file(s) ...');
      metadata.sort_values(by = options.key_field);

      with open(options.fasta_file, "r") as fasta_handle:
         data = SeqIO.parse(fasta_handle, "fasta");

         # Sort Fasta file so we organize upload, and can sync with metadata
         data = [f for f in sorted(data, key=lambda x : x.id)];

         # Splits into batches of options.batch (default 1000) or less records:
         splits = len(data)/options.batch;
         if splits < 1:
            splits = 1;
         data = np.array_split(data, splits);

         for count, sequences in enumerate(data):
            # Determine metadata rows pertinent to all sequences. They should be in same order
            id_index = [record.id for record in sequences];
            metabatch = metadata.loc[metadata[options.key_field].isin(id_index)];

            print('Files for ' + options.output_file + '.'+ str(count))
            # https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.to_csv.html
            if options.csv_file:
               metabatch.to_csv(options.output_file + '.'+ str(count) + '.id.csv', index=False); 
            else:
               # Or to tabular?  Issue of quoted strings?
               metabatch.to_csv(options.output_file + '.'+ str(count) + '.id.tsv', sep="\t", index=False);

            with open(options.output_file + '.'+ str(count) + '.id.fasta', 'w') as output_handle:
               SeqIO.write(sequences, output_handle, "fasta");



   """ 
   STEP 2: SUBMIT TO API

   Currently only virusseq () API is an option.
   API Information: https://github.com/cancogen-virus-seq/docs/wiki/How-to-Submit-Data-(API)

   Log in here to get API Key, good for 24 hours:
   https://portal.dev.cancogen.cancercollaboratory.org/login

   NOTE: Study_id field if not set to account project like "UHTC-ON" will
   trigger validation error "UNAUTHORIZED_FOR_STUDY_UPLOAD".
   """

   if options.api:

      if not options.api_token:
         sys.exit("An API user token is required for use with the [" + options.api + "] API.");

      # Retrieve batch files that need uploading
      batches = glob.glob("./" + options.output_file + '.*.id.fasta');

      if len(batches) > 0:
         print ('Performing batch upload ...');

      ################################### VirusSeq API ##########################
      # See: https://github.com/cancogen-virus-seq/docs/wiki/How-to-Submit-Data-(API)
      if options.api == 'VirusSeq_Portal':
         if options.development:  # TEST API ENDPOINT
            url = "https://muse.dev.cancogen.cancercollaboratory.org/";
         else:  # LIVE API ENDPOINT
            url = "https://muse.virusseq-dataportal.ca/";  

         custom_header = {'Authorization': 'Bearer ' + options.api_token}

         # TESTING: create an empty or junky .fasta and accompanying .tsv file
         # batches = ['test.0.id.fasta'];

         for filename in batches:
            filename_tsv = filename.replace('.id.fasta','.id.tsv');
            upload_files = [
               ('files', open(filename, 'rb')), 
               ('files', open(filename_tsv, 'rb'))
            ];
            print('Uploading batch: ' + filename);
            try:
               request = requests.post(url + 'submissions', files = upload_files, headers = custom_header);
            except Exception as err:
               sys.exit("API Server problem (check API URL?): " + repr(err));

            if request.status_code == 200:
               result = request.json();
               if ('submissionId' in result):
                  submission_id = result['submissionId'];
                  print('Batch was submitted! submissionId: ' + submission_id);

                  os.rename(filename, filename.replace('.id.','.'+submission_id+'.'));
                  os.rename(filename_tsv, filename_tsv.replace('.id.','.'+submission_id+'.'));

                  continue;    
               else:
                  print(result);
                  sys.exit("Resolve reported error, then rerun command!");
               

            # "Unauthorized client error status response" code occurs when key is not valid.
            # This halts processing of all remaining batches.
            if request.status_code == 401:
               print("Unauthorized client error status 401 response");
               sys.exit("Check to make sure your API key is current.");

            if request.status_code == 404:
               sys.exit("API service endpoint not recognized. Check API URL:" + url)
            
            request_error = request.json();
            status = request_error['status'];
            message = request_error['message'];

            print(request_error);
            errorInfo = request_error['errorInfo'];

            # "Bad Request" response status indicates something wrong with the input files.
            # We have json at this point.
            if request.status_code == 400: # or status == "BAD_REQUEST"
               print("Bad Request status 400 response.");

               if (status == "BAD_REQUEST"):

                  if message == 'Headers are incorrect!':
                     print (message);
                     print ("Unknown Headers:", errorInfo['unknownHeaders']);
                     print ("Missing Headers:", errorInfo['missingHeaders']);
                     sys.exit("Check to make sure the .tsv file headers are current.");

                  """
                  Example error:
                  "errorInfo":{"invalidFields":[{"fieldName":"specimen collector sample ID","value":"","reason":"NOT_ALLOWED_TO_BE_EMPTY","index":1},{"fieldName":"fasta header name","value":"","reason":"NOT_ALLOWED_TO_BE_EMPTY","index":1},{"fieldName":"study_id","value":" 23434","reason":"UNAUTHORIZED_FOR_STUDY_UPLOAD","index":1}]}}
                  """
                  if message == 'Found records with invalid fields':
                     for record in errorInfo['invalidFields']:
                        print ("row " + str(record['index']), '"' + record['fieldName'] + '"', 
                           record['reason'],
                           "value:",   record['value']
                        );
                     continue;

               # not sure where this should be positioned.
               # {"status":"FORBIDDEN","message":"Denied","errorInfo":{}}
               if (status == "FORBIDDEN"):
                  sys.exit("Your account associated with the API key has not been authorized, so this service is not available to you yet.");

            # Internal Server Error (code generated etc.)
            if request.status_code == 500:
               print(status);
               print(message);
               if message == "Flux#last() didn't observe any onNext signal":
                  print ("Does .tsv file have no data rows?");
               continue;

            print('Error: Unable to complete batch because of status code ' + str(request.status_code) + '\n' + request.text);
            continue;


   # STEP 3: Report on progress of each batch job that has been submitted.

      # Get list of batch files to get status for
      batches = glob.glob('./' + options.output_file + '.*.*.fasta');

      for filename in batches:
         
         submission_id = filename.split('.')[-2];

         if not submission_id == 'id':
            print ();
            print ('STATUS for: ' + filename);
            if options.short:
               error_max = options.short;
            else:
               error_max = options.batch;

            if (options.api == 'VirusSeq_Portal'):

               query = '?page=0&size=' + str(error_max) + '&sortDirection=ASC&sortField=submitterSampleId&submissionId=' + submission_id;

               if options.development:  # TEST API ENDPOINT
                  url = "https://muse.dev.cancogen.cancercollaboratory.org/";
               else:  # LIVE API ENDPOINT
                  url = "https://muse.virusseq-dataportal.ca/";  

               feedback = requests.get(url + 'uploads' + query, headers = custom_header);
               if feedback.status_code == 200:
                  response = feedback.json();

                  for submission in response['data']:
                      # print(submission);
                     if (submission['status'] == 'QUEUED'):
                        print (submission['submitterSampleId'], "Queued");

                     if (submission['status'] == 'ERROR'):
                        error_list = submission['error'].split('#')[1:];

                        old_item_label = '';

                        item_report = submission['submitterSampleId'];
                        for ptr, item in enumerate(error_list):
                           # Just show field name, not section
                           binding = item.split(':',1);
                           item_label = binding[0].split('/')[-1];
                           item_error = binding[1];
                           if (item_label == old_item_label):
                              item_report += item_error;
                           else:
                              item_report += '\n\t' + item_label + item_error;

                           old_item_label = item_label;

                        print (item_report);

               else:
                  print ('Status unavailable');

               print();


def log(log_handler, text): 
   print (text);
   log_handler.write(text);
   return text;


def get_metadata(log_handle, options):

   if not options.fasta_file:
      sys.exit(log(log_handle, "A sample sequencing fasta file is required."));

   if options.csv_file:
      metadata = pd.read_csv(options.csv_file, encoding = 'unicode_escape');
      file_suffix = 'csv';

   elif options.tsv_file:
      metadata = pd.read_table(options.tsv_file, delimiter='\t', encoding = 'unicode_escape');
      file_suffix = 'tsv';

   else:
      sys.exit(log(log_handle, "A sample contextual .tsv or .csv file is required."));

   # Determine match of .fasta file record's identifier field to sample metadata file column
   #if not options.key_field:
   #   options.key_field = metadata.columns[0];     # Defaults to first column of metadata

   if not options.key_field in metadata.columns:
      log(log_handle, 'The key field column you provided [' + options.key_field + '] was not found in the contextual data file\'s list of columns:');
      log(log_handle, metadata.columns);
      sys.exit(1);

   return metadata;
