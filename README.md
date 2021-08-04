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
- **The `--strict-size` `-s` flag should be used for batches with small numbers of samples**  
- **The `--max-full-batches` `-f` flag should be used to force as many batches as possible to have *exactly* `--max-n-per-batch` samples**  
- `--vars-to-balance` `-v` defines the list of variables for which more than one group should be present in each batch (default: `c('codedsiteid','randomgroupcode','sex_psca','older_than_40')`)   
- `--tissue-subset [TISSUE_CODE]` restricts balancing to a single tissue specified by `[TISSUE_CODE]`, which must match a value in either the 'Sample Type' column of one of `--shipment-manifest-excel` or the 'SampleTypeCode' column of one of `--api-metadata-csv`  
- `--block-randomization` can be used to include an "injection_order" column based on block randomization (i.e. individuals within a batch are randomized, keeping all samples from an individual together, and then samples within an indiviual are randomized)  
- `--separate-batch-files` can be used to write out one blinded batch assignment file *per batch* per tissue per assay instead of the default behavior of one file per tissue per assay  
- `--outdir` `-o` can be used to specify an output directory other than the current working directory  
- `--quietly` `-q` can be used to silence progress messages   

Expert level: 
- If you want the script to check more random combinations of samples before compromising the ideal batch sizes, increase `--max-inner-loop-iter` (default: 1e6)  
- If you want the script to check more random combinations of samples before decreasing the stringency of the batch balance checks, increase `--max-outer-loop-iter` (default: 1000 for `--strict-size` and 5000 otherwise)  
- `--balance-strictness` `-b` can be used to specify the initial strictness of the balance checks, with 10 being the strictest and 1 being the most lenient. By default, `--balance-strictness` is 1 when `--strict-size` is used and 10 otherwise  


## Outputs  
### Files  
By default, two files are written for each assay & tissue combination:  
- Blinded batch assignments in the format `files/precovid_[SAMPLE_TYPE]-samples_BLINDED-batch-assignments.csv` (see [example](examples/precovid_4-samples_BLINDED-batch-assignments.csv))  
- Unblinded batching metadata in the format `files/precovid_[SAMPLE_TYPE]-samples_UNBLINDED-batch-characteristics.csv`  

Use the `--separate-batch-files` flag to output separate blinded batch assignment files per batch, e.g. `files/precovid_[SAMPLE_TYPE]-samples_BLINDED-batch_3-assignments.csv`.  
 
### Plots 
One plot is saved for each assay & tissue combination. This plot includes the number of individuals and samples per batch as well as the balance across each level of each `--vars-to-balance`. **These plots should be visually examined to confirm that batches are adequately balanced, i.e. that numbers are reasonably distributed across each ROW.** (see [examples](examples/plots)).   

## Usage 

### Required R packages
```txt
data.table
readxl
testit
argparse
ggplot2
gtsummary
pheatmap
```

### Example commands 
Here is an example of how to run the script from the command line, assuming the shipment manifest Excel files and API metadata CSV files are in the same directory as this script. Include manifests and metadata for *all* pre-COVID clinical samples, i.e. both adult and pediatric shipments.  
```bash
Rscript randomization.R \
    --shipment-manifest-excel \
      Stanford_ADU830-10060_120720.xlsx \
      Stanford_PED830-10062_120720.xlsx \
    --api-metadata-csv \
      ADU830-10060.csv \
      PED830-10062.csv \
    --max-n-per-batch 94 \
    --outdir ~/Desktop/stanford_batches 
```  
Equivalently:  

```bash
Rscript randomization.R \
    -ship \
      Stanford_ADU830-10060_120720.xlsx \
      Stanford_PED830-10062_120720.xlsx \
    -api \
      ADU830-10060.csv \
      PED830-10062.csv \
    -max 94 \
    -o ~/Desktop/stanford_batches 
``` 
Remember to add the `--strict-size` or `-s` flag if the maximum number of samples per batch is small, e.g.:  
```bash
Rscript randomization.R \
    -ship ShipmentContents_BroadCarr_012521.xlsx \
    -api ADU822-10074.csv \
    -max 15 \
    -s \
    -o ~/Desktop/broad_batches 
```
Add the `--max-full-batches` or `-f` flag to force as many batches as possible to have *exactly* `--max-n-per-batch` samples, e.g.:  
```bash
Rscript randomization.R \
    -ship Stanford_ADU830-10060_120720.xlsx \
    -api ADU830-10060.csv \
    -max 88 \
    -o ~/Desktop/stanford_batches \
    --max-full-batches
```
To run the randomization script for a single tissue, use the `--tissue-subset` argument, where the supplied value must be a value in the 'Sample Type' column of one `--shipment-manifest-excel` OR a value in the 'SampleTypeCode' column of one `--api-metadata-csv`, e.g.:
```bash
Rscript randomization.R \
    -ship \
        Stanford_ADU830-10060_120720.xlsx \
        Stanford_PED830-10062_120720.xlsx \
    -api \
        ADU830-10060.csv \
        PED830-10062.csv \
    -max 94 \
    -o ~/Desktop/stanford_batches \
    --tissue-subset 06 \
    --overwrite
```
The `--overwrite` flag ignores existing batching outputs and overwrites the files. Without this flag, batching for a sample type will be skipped if a batching output already exists.  

To add an "injection_order" column based on block randomization, add the `--block-randomization` flag; to output separate blinded batch assignment files for each batch, add the `--separate-batch-files` flag:  
```bash
Rscript randomization.R \
    -ship \
        Stanford_ADU830-10060_120720.xlsx \
        Stanford_PED830-10062_120720.xlsx \
    -api \
        ADU830-10060.csv \
        PED830-10062.csv \
    -max 94 \
    -o ~/Desktop/stanford_batches \
    --tissue-subset 06 \
    --overwrite \
    --block-randomization \
    --separate-batch-files
```

See examples of the stdout for [large batches](examples/large-batches.out.log) and [small batches (`--strict-size`)](examples/small-batches.out.log).   

Alternatively, run the script interactively in RStudio by commenting out lines 19-66 and manually defining arguments below (see examples on lines 68-123), though this is not recommended.  

## Argument documentation
Run `Rscript randomization.R -h` to see this help message:  
```bash
usage: randomization.R [-h] -ship SHIPMENT_MANIFEST_EXCEL
                       [SHIPMENT_MANIFEST_EXCEL ...] -api API_METADATA_CSV
                       [API_METADATA_CSV ...] -max MAX_N_PER_BATCH [-s] [-f]
                       [-v VARS_TO_BALANCE] [-o OUTDIR] [-q]
                       [-inner MAX_INNER_LOOP_ITER]
                       [-outer MAX_OUTER_LOOP_ITER] [-b BALANCE_STRICTNESS]
                       [--overwrite] [--tissue-subset TISSUE_SUBSET]
                       [--block-randomization] [--separate-batch-files]

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
  -s, --strict-size     Force *all* batches to be as close to --max-n-per-
                        batch as possible. Most applicable for small batches
                        (e.g. < 20)
  -f, --max-full-batches
                        Force as many batches as possible to have *exactly*
                        --max-n-per-batch samples.
  -v VARS_TO_BALANCE, --vars-to-balance VARS_TO_BALANCE
                        Force batches to include samples from at least two
                        groups of each of these variables. Must be defined in
                        --api-metadata-csv
  -o OUTDIR, --outdir OUTDIR
                        Path to output directory
  -q, --quietly         Silence progress messages
  -inner MAX_INNER_LOOP_ITER, --max-inner-loop-iter MAX_INNER_LOOP_ITER
                        Max number of failed attempts to fit all samples in
                        batches before increasing the number of batches
  -outer MAX_OUTER_LOOP_ITER, --max-outer-loop-iter MAX_OUTER_LOOP_ITER
                        Max number of failed attempts to find optimally
                        balanced bacthes before relaxing the stringency of the
                        balance checks
  -b BALANCE_STRICTNESS, --balance-strictness BALANCE_STRICTNESS
                        Initial strictness of balance checks, with 10 being
                        the strictest and 1 being the most lenient
  --overwrite           Overwrite existing batching results
  --tissue-subset TISSUE_SUBSET
                        Run batching for a single tissue. Must be a value in
                        the 'Sample Type' column of one --shipment-manifest-
                        excel OR a value in the 'SampleTypeCode' column of one
                        --api-metadata-csv
  --block-randomization
                        Block randomization for metabolomics sites: samples
                        within a batch are ordered by individual; samples
                        within an individual are randomized. This adds an
                        'injection_order' column.
  --separate-batch-files
                        Write separate BLINDED output files per batch
```

## Help
For questions about the documentation or any issues with the code, please submit an issue or contact Nicole at nicolerg@stanford.edu. 
