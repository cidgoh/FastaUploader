""" fasta_batch_submit.py

Given a fasta file and a sample metadata file with a column that matches to
fasta file record identifiers, break both into respective sets of smaller
batches of records which are submitted to an API for processing.

Processing is two step: 

1) construct batches of files. Since two files are read and parsed in one go,
processing of them is reliable after that point, so no further error reporting
required after parsing.
 - Importantly, if rerunning, this step will be skipped unless -f --force 
 parameter is run.  But at moment input files are still required.

2) IF API option is included, submit each batch to API, wait for it to finish
or error out (capture error report) and proceed to next batch.

Authors: Damion Dooley, Nithu Sara John
Centre for Infectious Disease Epidemiology and One Health
August 24, 2021

Requires Biopython and Requests modules
- "pip install biopython" or https://biopython.org/wiki/Download
- "pip install requests"

Usage:
python fasta_batch_submit.py -c "20210713_AB_final set 1.csv" -f "consensus_renamed_final.fasta" -k "specimen collector sample ID"
"""

from Bio import SeqIO
import numpy as np
import optparse
import pandas as pd
import sys
import requests
import os, glob

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
   help="provide an output file name/path", default='output_');
parser.add_option('-k', '--key', dest="key_field",
   help="provide the metadata field name to match to fasta record identifier");
parser.add_option('-a', '--api', dest="api", default='virusseq',
   help="provide the target API to send data too.  A batch submission job will be initiated for it.");
parser.add_option('-u', '--user', dest="api_token",
   help="A user token is required for API access");
parser.add_option('-r', '--reset', dest="reset",
   help="Regenerate all batch files and clear all logs even if they already exist.");

options, args = parser.parse_args();

if not options.fasta_file:
   sys.exit("A sample sequencing fasta file is required.");

if options.csv_file:
   metadata = pd.read_csv(options.csv_file, encoding = 'unicode_escape');
   file_suffix = 'csv';

elif options.tsv_file:
   metadata = pd.read_table(options.tsv_file, delimiter='\t', encoding = 'unicode_escape');
   file_suffix = 'tsv';

else:
   sys.exit("A sample contextual .tsv or .csv file is required.");

# Determine match of .fasta file record's identifier field to sample metadata file column
#if not options.key_field:
#   options.key_field = metadata.columns[0];     # Defaults to first column of metadata

if not options.key_field in metadata.columns:
   print('The key field column you provided [' + options.key_field + '] was not found in the contextual data file\'s list of columns:')
   print(metadata.columns);
   sys.exit(1);

# If force, clear out all existing batch files and redo
if options.reset:
   # Add syntax check / security on options.output_file references?
   for filename in glob.glob("./" + options.output_file + '*'):
      os.remove(filename);

batches = glob.glob("./" + options.output_file + '*.fasta');

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

      # Splits into batches of 1000 or less records:
      data = np.array_split(data, len(data)/options.batch);

      for count, sequences in enumerate(data):
         # Determine metadata rows pertinent to all sequences. They should be in same order
         id_index = [record.id for record in sequences];
         metabatch = metadata.loc[metadata[options.key_field].isin(id_index)];

         print('Files for ' + options.output_file + str(count))
         # https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.to_csv.html
         if options.csv_file:
            metabatch.to_csv(options.output_file + str(count) + '.csv', index=False); 
         else:
            metabatch.to_csv(options.output_file + str(count) + '.tsv', sep="\t", index=False);

         with open(options.output_file + str(count) + '.fasta', 'w') as output_handle:
            SeqIO.write(sequences, output_handle, "fasta");

   batches = glob.glob("./" + options.output_file + '*.fasta');

""" Currently only virusseq () API is an option.
API Information: https://github.com/cancogen-virus-seq/docs/wiki/How-to-Submit-Data-(API)

Log in here to get API Key, good for 24 hours:
https://portal.dev.cancogen.cancercollaboratory.org/login

curl --location --request POST 'https://muse.virusseq-dataportal.ca/submissions' \
--header 'Authorization: Bearer eyJhbGciOiJSUzI1NiJ9.eyJpYXQiOjE2Mjk5MzE5ODQsImV4cCI6MTYyOTk0Mjc4NCwic3ViIjoiYTAwNDMwNjMtODJlNy00NDE2LWJjOGMtNGYyMjFjOTFjNGJiIiwiaXNzIjoiZWdvIiwianRpIjoiYzZmMTc3M2EtYjg2Yy00YWIzLTljNTYtZTIyYTRlYzcyZjdlIiwiY29udGV4dCI6eyJzY29wZSI6WyJNVVNFLVVIVEMtT04uUkVBRCIsIk1VU0UtVUhUQy1PTi5XUklURSJdLCJ1c2VyIjp7ImVtYWlsIjoiZGFtaW9uX2Rvb2xleUBzZnUuY2EiLCJzdGF0dXMiOiJBUFBST1ZFRCIsImZpcnN0TmFtZSI6IkRhbWlvbiIsImxhc3ROYW1lIjoiRG9vbGV5IiwiY3JlYXRlZEF0IjoxNjI5OTEwMDQ3Mjc0LCJsYXN0TG9naW4iOjE2Mjk5MzE5ODQ0NjQsInByZWZlcnJlZExhbmd1YWdlIjpudWxsLCJwcm92aWRlclR5cGUiOiJLRVlDTE9BSyIsInByb3ZpZGVyU3ViamVjdElkIjoiNDhhMTExNDMtOTAzNS00NTdjLTllYzgtNzA3MzZkZDNmYjU0IiwidHlwZSI6IlVTRVIiLCJncm91cHMiOlsiNGZkMTc2ZjAtMmE1NS00M2Q3LWJkMzItMWZhNmQwZWI4MTA5Il19fSwiYXVkIjpbXX0.kBqctB7Y3SibLihrzXA_A714HhZq2_pr_pYZ5EnF_QQQRy1_RfJrvxUylNrWcGZi-AWiqqW0WiFT9x8ypA481fZkP6QLikZbNDMrEcmXhtGvU0selqcFJqil-bH7GUTal5lLZdxAB6Ly5R05YcQnhscYrcZNs9teP86CIHbzpZ2pKhq27lKa5zDlBojdWcwnad6ANjHywKiW9vtDqPQxErN0A4TImxuP5PLk--aJyDvg22jq179uY4jUIrpg0Y0Y_jZJTiTlVbZcAUGE4moa3n7SqdeDw596mdcdEpmpJnRnW8A-AsISr2ghOWEXuMp6nPkyv5NcITBCbQgsDVOZbQ' \
--form 'files=@"test.fasta"' \
--form 'files=@"test.tsv"'

TESTING ON output_0.fasta and output_0.csv, knowing that it should be .TSV.

Authorized for MUSE-UHTC-ON
NOTE: Study_id field if not set to account project like "MUSE-UHTC-ON" will
trigger validation error "UNAUTHORIZED_FOR_STUDY_UPLOAD".

"""

if options.api:
   print ('Performing batch upload ...');
   if not options.api_token:
      sys.exit("An API user token is required for use with the [" + options.api + "] API.");

   ################################### VirusSeq API ##########################
   if options.api == 'virusseq':
      url = "https://muse.virusseq-dataportal.ca/submissions";

      # https://www.kite.com/python/answers/how-to-add-custom-headers-to-a-post-request-in-python
      custom_header = {'Authorization': 'Bearer ' + options.api_token}

      #TESTING
      batches = ['test.fasta'];

      for filename in batches:

         upload_files = [
            ('files', open(filename, 'rb')), 
            ('files', open(filename.replace('.fasta','.tsv'), 'rb'))
         ];
         print('Processing batch: ' + filename);
         request = requests.post(url, files = upload_files, headers = custom_header);

         if request.status_code == 200:
            print('Processed?');
            result = request.json();
            print(result);
            continue            

         err = request.text;

         # "Unauthorized client error status response" code occurs when key is not valid.
         # This halts processing of all remaining batches.
         if request.status_code == 401:
            print("Unauthorized client error status 401 response");
            sys.exit("Check to make sure your API key is current.");

         request_error = request.json();
         status = request_error['status'];
         message = request_error['message'];
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


               if message == 'Found records with invalid fields':
                  for record in errorInfo['invalidFields']:
                     print ("row " + str(record['index']), '"' + record['fieldName'] + '"', 
                        record['reason'],
                        "value:",   record['value']
                     );
                  continue;

               """
               "errorInfo":{"invalidFields":[{"fieldName":"specimen collector sample ID","value":"","reason":"NOT_ALLOWED_TO_BE_EMPTY","index":1},{"fieldName":"fasta header name","value":"","reason":"NOT_ALLOWED_TO_BE_EMPTY","index":1},{"fieldName":"study_id","value":" 23434","reason":"UNAUTHORIZED_FOR_STUDY_UPLOAD","index":1}]}}
               """

            # not sure where this is positioned.
            # {"status":"FORBIDDEN","message":"Denied","errorInfo":{}}
            if (status == "FORBIDDEN"):
               sys.exit("Your account has not been authorized for use yet for file upload.");

         # Internal Server Error (code generated etc.)
         if request.status_code == 500:
            print(status);
            print(message);
            if message == "Flux#last() didn't observe any onNext signal":
               print ("Does .tsv file have no data rows?");
            continue;

         print('Error: Unable to complete batch because of status code ' + str(request.status_code) + '\n' + request.text);


         # except Exception as e: 
         #   print ("No content provided in response from\n" + url + "\n Is server down?");
         #   print(e);
