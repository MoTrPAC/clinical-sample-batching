# MoTrPAC clinical sample batching 
Contact: Nicole Gay (nicolerg@stanford.edu) 

Use [`randomization.R`](randomization.R) to make well-balanced batches of MoTrPAC human samples in terms of clinical site, intervention group, age, and sex. 

## Required inputs  
- Shipment manifest Excel file(s) from the Biorepository, e.g. `Stanford_ADU830-10060_120720.xlsx` (`--shipment-manifest-excel` `-ship`)   
- Corresponding CSV file(s) from the Biospecimen Metadata Download API on [motrpac.org](https://www.motrpac.org/), e.g. `ADU830-10060.csv` (`--api-metadata-csv` `-api`)   
- Maximum number of samples per batch (`--max-n-per-batch` `-max`)  

> IMPORTANT: Include manifests for both adult *and* pediatric samples to randomize studies together.  

> IMPORTANT: If the manifests include multiple aliquots of the same samples for different assays, you **MUST** add an `assay` column to *either* of these input files to distiguish the different sets of aliquots. For example, if muscle samples are being processed for both ATAC-seq and RNA-seq at Stanford, add an `assay` column with the values `rnaseq` and `atacseq`. The values themselves do not matter as long as they separate the sets of aliquots.  

## Optional arguments  
- **The `--strict-size` `-t` flag should be used for batches with small numbers of samples.**  
- `--vars-to-balance` `-v` defines the list of variables for which more than one group should be present in each batch (default: `c('codedsiteid','randomgroupcode','sex_psca','older_than_40')`). 
- `--outdir` `-o` can be used to specify an output directory other than the current working directory.  
- `--quietly` `-q` can be used to silence progress messages.  

## Outputs  
- Two files per assay & tissue combination:  
  - Blinded batch assignments in the format `precovid_[SAMPLE_TYPE]-samples_BLINDED-batch-assignments.csv` (see example [here](examples/precovid_4-samples_BLINDED-batch-assignments.csv))  
  - Unblinded batching metadata in the format `precovid_[SAMPLE_TYPE]-samples_UNBLINDED-batch-characteristics.csv`  
- Summary of batch characteristics (`stdout`) (see example [here](examples/out.log))  

## Usage 

### Required R packages
```txt
data.table
readxl
testit
argparse
```

### Example commands 
Here is an example of how to run the script from the command line, assuming the shipment manifest Excel files and API metadata CSV files are in the same directory as this script. Include manifests and metadata for *all* pre-COVID clinical samples, i.e. both adult and pediatric shipments.  
```bash
Rscript randomization.R \
    --shipment-manifest-excel Stanford_ADU830-10060_120720.xlsx Stanford_PED830-10062_120720.xlsx \
    --api-metadata-csv ADU830-10060.csv PED830-10062.csv \
    --max-n-per-batch 94 \
    --outdir ../batches 
```  
Equivalently:  

```bash
Rscript randomization.R \
    -ship Stanford_ADU830-10060_120720.xlsx Stanford_PED830-10062_120720.xlsx \
    -api ADU830-10060.csv PED830-10062.csv \
    -max 94 \
    -o ../batches 
```  
A summary of batching statistics is printed to the console. To save all output to a log file for later reference, add ` > out.log 2>&1` to the end of the command, e.g.: 
```bash
Rscript randomization.R \
    -ship Stanford_ADU830-10060_120720.xlsx Stanford_PED830-10062_120720.xlsx \
    -api ADU830-10060.csv PED830-10062.csv \
    -max 94 \
    -o ../batches > ../batches/out.log 2>&1
```
See an example of this log file [here](examples/out.log). 

Alternatively, run the script interactively in RStudio by commenting out lines 15-37 and manually defining arguments below (see examples on lines 39-55).  

## Argument documentation
Run `Rscript randomization.R -h` to see this help message:  
```bash
usage: randomization.R [-h] -ship SHIPMENT_MANIFEST_EXCEL
                       [SHIPMENT_MANIFEST_EXCEL ...] -api API_METADATA_CSV
                       [API_METADATA_CSV ...] -max MAX_N_PER_BATCH [-s]
                       [-v VARS_TO_BALANCE] [-o OUTDIR] [-q]

optional arguments:
  -h, --help            show this help message and exit
  -ship SHIPMENT_MANIFEST_EXCEL [SHIPMENT_MANIFEST_EXCEL ...], --shipment-manifest-excel SHIPMENT_MANIFEST_EXCEL [SHIPMENT_MANIFEST_EXCEL ...]
                        Path(s) to shipment manifest Excel files, e.g.
                        Stanford_ADU830-10060_120720.xlsx
                        Stanford_PED830-10062_120720.xlsx
  -api API_METADATA_CSV [API_METADATA_CSV ...], --api-metadata-csv API_METADATA_CSV [API_METADATA_CSV ...]
                        Path(s) to sample metadata from web API, e.g.
                        ADU830-10060.csv PED830-10062.csv
  -max MAX_N_PER_BATCH, --max-n-per-batch MAX_N_PER_BATCH
                        Max number of samples per batch
  -s, --strict-size     Force all batches to be as close to --max-n-per-batch
                        as possible. Most applicable for small batches (e.g. <
                        20)
  -v VARS_TO_BALANCE, --vars-to-balance VARS_TO_BALANCE
                        Force batches to include samples from at least two
                        groups of each of these variables. Must be defined in
                        --api-metadata-csv
  -o OUTDIR, --outdir OUTDIR
                        Path to output directory
  -q, --quietly         Silence progress messages
```

## Help
For questions about the documentation or any issues with the code, please submit an issue or contact Nicole at nicolerg@stanford.edu. 
