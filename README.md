# fasta-harmonizer
A tool for batch processing large fasta files and accompanying metadata table to repositories via API

The python **fasta_batch_submit.py** script breaks large fasta files (e.g. 500mb) and related (one-to-one) tabular contextual data into smaller batches of 1000 or some specified # of records which can then be uploaded to a given sequence repository if an API endpoint is selected.  Currently API interface is developed for VirusSeq (https://virusseq-dataportal.ca/). 

This tool is developed by the SFU Centre for Infectious Disease Epidemiology and One Health in conjunction with VirusSeq and it works well with [DataHarmonizer](https://github.com/Public-Health-Bioinformatics/DataHarmonizer)!

Authors: Damion Dooley, Nithu Sara John

## Details 

Given a fasta file and a sample metadata file with a column that matches to fasta file record identifiers, break both into respective sets of smaller batches of records which are submitted to an API for processing.

Processing is two step: 

1) Construct batches of files. Since only two files are read and parsed in one go,
processing of them is reliable after that point, so no further error reporting
is needed during the batch file generation process.
   1) Importantly, if rerunning fasta_batch_submit.py, this step will be skipped unless -f --force parameter is run.  Currently input files are still required in this case.

1) IF API option is included, submit each batch to API, wait for it to finish
or error out (capture error report) and proceed to next batch. 
   1) Some types of error trigger sudden death, i.e. sys.exit() because they would apply to any subsequent API batch calls.  For example missing tabular data column names will trigger an exit. Once resolved, rerun with -f to force regeneration of output files.

Requires Biopython and Requests modules

- "pip install biopython"
- "pip install requests"

## Usage
python fasta_batch_submit.py [options]

Options:

    -h, --help
      show this help message and exit
    -f FASTA_FILE, --fasta=FASTA_FILE
      provide a fasta file name
    -c CSV_FILE, --csv=CSV_FILE
      provide a COMMA delimited sample contextual data file name
    -t TSV_FILE, --tsv=TSV_FILE
      provide a TAB delimited sample contextual data file name
    -b BATCH, --batch=BATCH
      provide number of fasta records to include in each batch
    -o OUTPUT_FILE, --output=OUTPUT_FILE
      provide an output file name/path
    -k KEY_FIELD, --key=KEY_FIELD
      provide the metadata field name to match to fasta record identifier
    -a API, --api=API     
      provide the target API to send data too.  A batch submission job will be initiated for it. Default is "VirusSeq_Portal"
    -u API_TOKEN, --user=API_TOKEN
      an API user token is required for API access
    -r, --reset
      Regenerate all batch files and begin API resubmission process even if batch files already exist under given output file pattern.
                        
For example:

python fasta_batch_submit.py -c "final set 1.csv" -f "consensus_final.fasta" -k "specimen collector sample ID" -u [enter API key here]


