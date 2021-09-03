# Fasta Uploader
A tool for batch processing large fasta files and accompanying metadata table to repositories via API

The python **fasta_uploader.py** script breaks large fasta files (e.g. 500mb) and related (one-to-one) tab-delimited sample contextual data into smaller batches of 1000 or some specified # of records which can then be uploaded to a given sequence repository if an API endpoint is selected.  Currently there is one option for the API interface: [VirusSeq](https://virusseq-dataportal.ca/). 

This tool is developed by the SFU Centre for Infectious Disease Epidemiology and One Health in conjunction with VirusSeq and it works well with [DataHarmonizer](https://github.com/Public-Health-Bioinformatics/DataHarmonizer)!

Authors: Damion Dooley, Nithu Sara John

## Details 

Given a fasta file and a sample metadata file with a column that matches to fasta file record identifiers, break both into respective sets of smaller batches of records which are submitted to an API for processing.

Processing is three step: 

1) Construct batches of files. Since only two files are read and parsed in one go,
processing of them is reliable after that point, so no further error reporting
is needed during the batch file generation process.
   1) Importantly, if rerunning fasta_batch_submit.py, this step will be skipped unless -f --force parameter is run.  Currently input files are still required in this case.

1) IF API option is included, submit each batch to API, wait for it to finish
or error out (capture error report) and proceed to next batch. 
   1) Some types of error trigger sudden death, i.e. sys.exit() because they would also occur in subsequent API batch calls.  For example missing tabular data column names will trigger an exit. Once resolved, rerun with -f to force regeneration of output files.

1) Report on processing status of existing API requests.  Some may be queued, others may have been processed successfully, and others may have line-by-line errors in field content that must be addressed.  Revise given batch files and run again.

Requires Biopython and Requests modules

- "pip install biopython"
- "pip install requests"

## Usage
Run the command in a folder with the appropriate input files, and output files can be generated there too.  Rerun it in the same folder to incrementally fix any submission errors and then restart submission.

**python fasta_uploader.py [options]**

### Options:

    -h, --help
      Show this help message and exit.
    -f FASTA_FILE, --fasta=FASTA_FILE
      Provide a fasta file name.
    -c CSV_FILE, --csv=CSV_FILE
      Provide a COMMA delimited sample contextual data file name.
    -t TSV_FILE, --tsv=TSV_FILE
      Provide a TAB delimited sample contextual data file name.
    -b BATCH, --batch=BATCH
      Provide number of fasta records to include in each batch. Default is 1000.
    -o OUTPUT_FILE, --output=OUTPUT_FILE
      Provide an output file name/path.
    -k KEY_FIELD, --key=KEY_FIELD
      Provide the metadata field name to match to fasta record identifier.
    -n BATCH_NUMBER, --number=BATCH_NUMBER
      Process only given batch number to API instead of all batches.
      
Parameters involved in optional API call:
      
    -a API, --api=API     
      Provide the target API to send data too.  A batch submission job will be initiated for it. Default is "VirusSeq_Portal".
    -u API_TOKEN, --user=API_TOKEN
      An API user token is required for API access.
    -d, --dev
      Test against a development server rather than live one.  Provide an API endpoint URL.
     -s, --short
      Report up to given # of fasta record related errors for each batch submission.  Useful for taking care of repeated errors first based on first instance.
     -r, --reset
      Regenerate all batch files and begin API resubmission process even if batch files already exist under given output file pattern.
   
For example:

python fasta_uploader.py -f "consensus_final.fasta" -c "final set 1.csv" -k "fasta header name" -u [enter API key here] -s 2


