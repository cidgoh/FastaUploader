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

1) IF API option is included, submit each *.queued.fasta batch to API, wait for it to finish
or error out (capture error report) and proceed to next batch. 
   1) There is an option to just try submitting one of the batches, e.g. the first one, via "-n 0" parameter.  This allows error debugging of just the first batch.  Once error patterns are determined, those that apply to remaining source contextual data can be applied, and first batch removed from that file, and the whole batching can be redone using -r reset.
   1) Some types of error trigger sudden death, i.e. sys.exit() because they would also occur in subsequent API batch calls.  For example missing tabular data column names will trigger an exit. Once resolved, rerun with -f to force regeneration of output files.

1) Report on processing status of existing API requests.  Some may be queued, others may have been processed successfully, and others may have line-by-line errors in field content that are turned into new [output file batch.#].queued.fasta and [output batch.#].queued.tsv files which can be edited and then submitted by rerunning the program.

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
    -m METADATA_FILE, --metadata=METADATA_FILE
      Provide a COMMA .csv or TAB .tsv delimited sample contextual data file name.
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

    python fasta_uploader.py -f "consensus_final.fasta" -m "final set 1.csv" -k "fasta header name" -a VirusSeq_Portal -u ENTER_API_KEY_HERE

This will convert consensus_final.fasta and related final set 1.csv contextual data records into batches of 1000 records by default, and will begin submitting each batch to the VirusSeq portal.

    python fasta_uploader.py -f "consensus_final.fasta" -m "final set 1.csv" -n 0 -k "fasta header name" -a VirusSeq_Portal -u ENTER_API_KEY_HERE

Like the above but only first batch is submitted so that one can see any errors, and if they apply to all batches, can fix them in original "final set 1.csv" file. Once batch 0 is fixed, all its records can be removed from the consensus_final.fasta and final set 1.csv.csv files, and the whole job can be resubmitted.
