#' Randomize samples into balanced batches
#' 
#' Given different batching considerations for all MoTrPAC Chemical Analysis
#' Sites (CAS), this function generates batches of a given size that are 
#' balanced in terms of the specified phenotypic variables. Outputs are written
#' to a local directory. 
#'
#' @param shipment_manifest character vector, path(s) or URL to shipment manifests, 
#'   e.g. \code{c('Stanford_ADU830-10060_120720.xlsx', 'Stanford_PED830-10062_120720.xlsx')}.
#'   Can be Excel files (.xlsx), tab-delimited text (.tsv), or comma-separated text (.csv).  
#'   Must include a "viallabel" column (case insensitive) to map to \code{api_metadata_csv}. 
#'   See details for other expected columns. 
#' @param api_metadata character vector, path(s) or URL to sample metadata from web API, 
#'   e.g. \code{'ADU830-10060.csv', 'PED830-10062.csv'}.
#'   Can be Excel files (.xlsx), tab-delimited text (.tsv), or comma-separated text (.csv).  
#'   Must include a "viallabel" column (case insensitive) to map to \code{shipment_manifest_excel}. 
#'   See details for other expected columns. 
#' @param max_n_per_batch integer, max number of samples per batch
#' @param strict_size bool, force *all* batches to be as close to \code{max_n_per_batch} as possible. 
#'   Most applicable for small batches (e.g. < 20). Default: FALSE
#' @param max_full_batches bool, whether to force as many batches as possible to have *exactly* \code{max_n_per_batch} samples.
#'   Default: FALSE 
#' @param vars_to_balance character vector, force batches to include samples from at 
#'   least one level of each of these variables. Must be defined in the column names of 
#'   \code{api_metadata_csv} or \code{shipment_manifest_excel} (case insensitive). 
#'   Default: \code{c('codedsiteid','randomgroupcode','sex_psca','older_than_median')} 
#' @param outdir character, path to output directory. Current working directory by default. 
#' @param verbose bool, whether to print progress messages. Default: TRUE. 
#' @param max_inner_loop_iter integer, max number of failed attempts to fit all samples 
#'   in batches before increasing the number of batches. Decreasing this number 
#'   may generate batches more quickly, but the number of samples per batch may not be optimal.
#' @param max_outer_loop_iter integer, max number of failed attempts to find optimally 
#'   balanced batches before relaxing the stringency of the balance checks.
#'   Decreasing this number may generate batches more quickly, but batches may not
#'   be optimally balanced. 
#' @param overwrite bool, whether to overwrite existing batching results.
#'   Defaults to current working directory. 
#' @param balance_strictness integer, initial strictness of balance checks, 
#'   with 10 being the strictest and 1 being the most lenient. Defaults to 1 if \code{strict_size = TRUE}
#'   and 9 otherwise. Decreasing this number may generate batches more quickly, but batches may not
#'   be optimally balanced. 
#' @param tissue_subset character, only specified to run batching for a single tissue. 
#'   Must be a value in the 'Sample Type' column of one \code{shipment_manifest_excel} 
#'   OR a value in the 'SampleTypeCode' column of one \code{api_metadata_csv}."
#' @param block_randomization bool, whether to use block randomization. 
#'   Samples within a batch are ordered by individual; 
#'   samples within an individual are randomized. This adds an 'injection_order' column.
#'   Default: FALSE
#' @param separate_batch_files bool, whether to write separate BLINDED output files per batch.
#'   Default: FALSE 
#'
#' @export
#'
#' @examples
#' # Get paths to shipment and metadata files.
#' # Here we use examples provided in this R package. 
#' # These are examples of real MoTrPAC clinical sample batching files, where 
#' # PIDs, vial labels, and barcodes were changed to random numbers. 
#' # Note that both adult and pediatric samples are included in these examples,
#' # which means pediatric and adult samples are randomized together. 
#' shipment = system.file("extdata", "shipment-manifest.xlsx", package = "SampleBatchR")
#' metadata = system.file("extdata", "biospecimen-metadata.csv", package = "SampleBatchR")
#' 
#' # For each tissue type present in the input files,
#' # generate batches with max 96 samples per batch.
#' # Save outputs to Desktop.
#' \dontrun{
#' batch_samples(
#'  shipment, 
#'  metadata, 
#'  max_n_per_batch = 96,
#'  outdir = "~/Desktop/batches/all" 
#' )
#' }
#' 
#' # Generate batches of max 96 samples for muscle samples only (code 6).
#' # Save outputs to Desktop.
#' \dontrun{
#' batch_samples(
#'  shipment, 
#'  metadata, 
#'  max_n_per_batch = 96,
#'  outdir = "~/Desktop/batches/muscle",
#'  tissue_subset = "6"
#' )
#' }
#' 
#' # Generate small batches of max 15 samples.
#' # For small batch sizes, `strict_size` should be set to TRUE, 
#' # which forces as many batches as possible to have exactly N samples. 
#' \dontrun{
#' batch_samples(
#'  shipment, 
#'  metadata, 
#'  max_n_per_batch = 15,
#'  outdir = "~/Desktop/batches/small",
#'  strict_size = TRUE, 
#'  tissue_subset = "6"
#' )
#' }
#' 
#' # To force as many batches as possible to have exactly N samples
#' # for larger batch sizes, e.g., >= 30, set `max_full_batches` to TRUE
#' \dontrun{
#' batch_samples(
#'  shipment, 
#'  metadata, 
#'  max_n_per_batch = 88,
#'  outdir = "~/Desktop/batches/max_full",
#'  max_full_batches = TRUE, 
#'  tissue_subset = "6"
#' )
#' }
#' 
#' # To add an "injection_order" column based on block randomization,
#' # set `block_randomization` to TRUE. For block randomization, samples within 
#' # a batch are ordered by individual; samples within an individual are randomized. 
#' # Set `overwrite` to TRUE to force overwrite of existing results. 
#' \dontrun{
#' batch_samples(
#'  shipment, 
#'  metadata, 
#'  max_n_per_batch = 88,
#'  outdir = "~/Desktop/batches/block_rand",
#'  block_randomization = TRUE, 
#'  tissue_subset = "6",
#'  overwrite = TRUE
#' )
#' }
#' 
#' @details 
#' This function is currently written for compatibility with MoTrPAC clinical sample
#' shipment manifests and biospecimen metadata files. It can be further generalized if
#' there is demand to use this function for other projects. 
#' 
#' **Inputs:**  
#' 
#' Expected column names in the input (case insensitive):  
#' * Required: \code{viallabel} in both \code{shipment_manifest} and \code{api_metadata}. Sample-level identifier. Used to merge these files. 
#' * Required: \code{barcode} in either \code{shipment_manifest} or \code{api_metadata}. Another sample-level identifier. Can be NA. 
#' * Required: \code{pid} in either \code{shipment_manifest} or \code{api_metadata}. Participant-level identifier. 
#' If there is only one sample per participant, this column can be a copy of \code{viallabel}. 
#' * Required: Each variable listed in \code{vars_to_balance} must be a column name in either \code{shipment_manifest} or \code{api_metadata}.
#' * Required: \code{box} in either \code{shipment_manifest} or \code{api_metadata}. Box identifier to indicate which 
#' box/plate the sample was shipped in. 
#' * Required: \code{position} in either \code{shipment_manifest} or \code{api_metadata}. Provides the position
#' of the sample in the box/plate in which it was shipped. 
#' * Required if \code{tissue_subet} is defined: \code{sampletypecode} in \code{api_metadata} or \code{sample type} in \code{shipment_manifest}.  
#' * Optional: \code{assay} in either \code{shipment_manifest} or \code{api_metadata}. Define this column if subsets of samples are destined for different assays.
#' * Required to generate a heatmap of batch balance: \code{'calculatedage','sex_psca','codedsiteid','randomgroupcode'}  
#' 
#' **Outputs:**  
#' 
#' This function writes outputs to a local directory, specified with the \code{outdir} argument.
#' 
#' *Files*  
#' 
#' By default, two files are written for each assay & tissue combination:
#' * Blinded batch assignments in the format \code{files/[SAMPLE_TYPE]-samples_BLINDED-batch-assignments.csv}  
#' * Unblinded batching metadata in the format \code{files/[SAMPLE_TYPE]-samples_UNBLINDED-batch-characteristics.csv}  
#' 
#' Set \code{separate_batch_files} to TRUE to output separate blinded batch assignment files per batch, 
#' e.g. \code{files/[SAMPLE_TYPE]-samples_BLINDED-batch_3-assignments.csv}.  
#' 
#' *Plots*  
#' 
#' One heatmap is saved for each assay & tissue combination. This plot includes the number of 
#' individuals and samples per batch as well as the balance across each level of each \code{vars_to_balance}. 
#' These plots should be visually examined to confirm that batches are adequately balanced, 
#' i.e. that numbers are reasonably distributed across each ROW. Note that these plots
#' are currently only generated if the following columns are included in the input files:
#' 'calculatedage','sex_psca','codedsiteid','randomgroupcode'. If these columns are not
#' included, then batch-by-level tables are printed to the console for each \code{vars_to_balance}. 
#' 
batch_samples = function(shipment_manifest,
                         api_metadata,
                         max_n_per_batch,
                         strict_size = FALSE,
                         max_full_batches = FALSE,
                         vars_to_balance = c('codedsiteid','randomgroupcode','sex_psca','older_than_median'),
                         outdir = ".",
                         verbose = TRUE, 
                         max_inner_loop_iter = 1e5,
                         max_outer_loop_iter = NULL,
                         overwrite = FALSE,
                         balance_strictness = NULL,
                         tissue_subset = NULL,
                         block_randomization = FALSE,
                         separate_batch_files = FALSE){
  
  shipments = shipment_manifest
  apis = api_metadata
  max_full = max_full_batches
  balance_vars = vars_to_balance
  init_balance_strictness = balance_strictness
  if(!is.null(tissue_subset)) tissue_subset = as.character(tissue_subset)
  
  # initialize balance strictness and enforce bounds 
  if(is.null(init_balance_strictness)){
    if(strict_size){
      init_balance_strictness = 1
    }else{
      init_balance_strictness = 9
    }
  }
  if(init_balance_strictness > 10){
    init_balance_strictness = 10
  }
  if(init_balance_strictness < 1){
    init_balance_strictness = 1
  }
  
  if(max_full & strict_size){
    stop(paste("This script is not currently designed for cases where *both* 'max_full_batches' and 'strict_size' are TRUE.",
               "Please set only one of the two options to TRUE. 'strict_size' should be used for small batches, e.g. 'max_n_per_batch' < 20;",
               "'max_full_batches' should be used for larger batches where the user wants to maximize the number of batches with exactly", 
               "'max_n_per_batch' samples."))
  }
  
  if(is.null(max_outer_loop_iter)){
    if(strict_size){
      max_outer_loop_iter = 1000
    }else{
      max_outer_loop_iter = 5000 
    }
  }
  
  # check formats
  patterns = c("\\.xlsx$","\\.csv$","\\.tsv$")
  for(f in c(shipments,apis)){
    if(!any(sapply(patterns, function(p) grepl(p, f, ignore.case=T)))){
      stop("At least one of the input files is not in one of the accepted formats: .xlsx, .csv, .tsv.")
    }
  }
  
  #### load manifests ####
  
  ship_list = list()
  for (s in shipments){
    ship_list[[s]] = read_input(s)
  }
  ship = rbindlist(ship_list, fill=T)
  
  # read in API metadata 
  api_list = list()
  for (a in apis){
    api_list[[a]] = read_input(a)
  }
  api = rbindlist(api_list, fill=T)
  
  # make colnames lowercase for simplicity
  colnames(ship) = tolower(colnames(ship))
  colnames(api) = tolower(colnames(api))
  balance_vars = tolower(balance_vars)
  
  # binarize age on median 
  med = median(api[,calculatedage])
  api[,older_than_median := calculatedage > med]
  
  # merge
  api[,viallabel := as.character(viallabel)]
  ship[,viallabel := as.character(viallabel)]
  #ship[,`2d barcode` := as.character(`2d barcode`)]
  # PBMCs don't have barcodes?
  all_meta = merge(api, ship, by='viallabel')
  stopifnot(nrow(all_meta) == nrow(ship))
  all_meta[all_meta=='.'] = NA
  
  if(nrow(all_meta)!=nrow(ship)){
    warning(sprintf("The number of samples in the merged metadata (%s) is not the same as the number of samples in the shipment manifest '%s' (%s). Check for a merging error.\n",
                    nrow(all_meta),
                    paste0(basename(shipments), collapse=', '),
                    nrow(ship)))
  }
  if(nrow(all_meta)!=nrow(api)){
    warning(sprintf("The number of samples in the merged metadata (%s) is not the same as the number of samples in the biospecimen metadata '%s' (%s). Check for a merging error.\n",
                    nrow(all_meta),
                    paste0(basename(apis), collapse=', '),
                    nrow(api)))
  }
  
  # check if we need to subset by a tissue
  if(!is.null(tissue_subset)){
    
    found_match = FALSE
    # try to find a matching value in either `sample type` or `sampletypecode`
    if("sample type" %in% colnames(all_meta) | "sampletypecode" %in% colnames(all_meta)){
      
      # add dummy columns 
      if(!"sample type" %in% colnames(all_meta)){
        all_meta[,`sample type` := NA_character_]
      }else if(!"sampletypecode" %in% colnames(all_meta)){
        all_meta[,sampletypecode := NA_character_]
      }
      # look for matches
      if(tissue_subset %in% all_meta[,sampletypecode]){
        all_meta = all_meta[sampletypecode == tissue_subset]
        found_match = TRUE
      }else if(tissue_subset %in% all_meta[,`sample type`]){
        all_meta = all_meta[`sample type` == tissue_subset]
        found_match = TRUE
      }else if(tolower(tissue_subset) %in% tolower(all_meta[,sampletypecode])){
        all_meta = all_meta[tolower(sampletypecode) == tolower(tissue_subset)]
        found_match = TRUE
      }else if(tolower(tissue_subset) %in% tolower(all_meta[,`sample type`])){
        all_meta = all_meta[tolower(`sample type`) == tolower(tissue_subset)]
        found_match = TRUE
      }else if(is.numeric(all_meta[,`sample type`])){
        if(as.numeric(tissue_subset) %in% all_meta[,`sample type`]){
          all_meta = all_meta[`sample type` == as.numeric(tissue_subset)]
          found_match = TRUE
        }
      }else if(is.numeric(all_meta[,sampletypecode])){
        if(as.numeric(tissue_subset) %in% all_meta[,sampletypecode]){
          all_meta = all_meta[sampletypecode == as.numeric(tissue_subset)]
          found_match = TRUE
        }
      }
    }else{
      stop("Neither 'sample type' nor 'sampletypecode' were found in column names of the input, but 'tissue_subset' is specified.")
    }
    
    if(!found_match){
      stop(sprintf("Unable to match tissue_subset '%s' to a value in either 'api_metadata$SampleTypeCode' or 'shipment_manifest$`Sample Type`'. Available options to subset by tissue: %s",
                   tissue_subset, 
                   paste0(c(unique(all_meta[,sampletypecode]), unique(all_meta[,`sample type`])), collapse=', ')))
    }
  }
  
  # make sure randomization variables exist
  missing = c()
  for(b in tolower(balance_vars)){
    if(!b %in% colnames(all_meta)){
      missing = c(missing, b)
    }
  }
  if(length(missing) > 0){
    warning(sprintf("The following variables were specified in 'vars_to_balance' but were not found in the input files:\n %s",
                    paste0(missing, collapse=", ")))
    balance_vars = balance_vars[!balance_vars %in% missing]
  }
  if(length(balance_vars) == 0){
    warning("No variables will be checked for balance across batches. Batching samples anyway.")
  }
  
  # subset to existing samples
  all_meta = all_meta[!(is.na(viallabel) | is.na(box))]
  
  #### setup ##### 
  
  # make outdirs
  if(!dir.exists(sprintf("%s/plots", outdir))){
    dir.create(sprintf("%s/plots", outdir), recursive = TRUE)
  }
  if(!dir.exists(sprintf("%s/files", outdir))){
    dir.create(sprintf("%s/files", outdir), recursive = TRUE)
  }
  
  # "study" is redundant with "randomgroupcode"
  # there can be multiple bid per pid. pid id human participant id. use pid as identifier
  # visitcode = baseline versus post. ignore this because all samples from an individual will be together
  
  # make column names more flexible 
  cols = c("viallabel", "barcode", "pid", "box", "position", balance_vars)
  for(var in c("assay","sampletypecode","codedsiteid","randomgroupcode","sex_psca","calculatedage","older_than_median")){
    if(var %in% colnames(all_meta)){
      cols = c(cols, var)
    }
  }
  cols = unique(cols)
  all_meta = all_meta[,cols, with=FALSE]
  
  # note columns with 0 variance
  if(verbose){
    novar = c()
    for(c in colnames(all_meta)){
      if(length(unique(all_meta[,get(c)]))==1){
        novar = c(novar, c)
      }
    }
    if(length(novar) > 0){
      message(sprintf("These columns have 0 variance in the merged metadata: %s",
                      paste0(novar, collapse=', ')))
    }
  }

  #table(all_meta[,sampletypecode], all_meta[,randomgroupcode])
  # 01 = Human Serum
  # 02 = Human EDTA Plasma
  # 03 = Human EDTA Packed Cells
  # 04 = Human PAXgene RNA 05 = Human PBMC
  # 06 = Human Muscle
  # 07 = Human Adipose
  # 08 = Human Urine
  # 09 = Human Stool
  # 10 = Human Muscle Powder
  # 11 = Human Adipose Powder
  # 12 = Human EDTA 2 Plasma
  # 13 = Human Heparin Plasma
  # 14 = Human IMAT
  # 15 = Human Packed Cells DMSO 
  # 16 = Human Muscle Cell Passage 2 
  # 17 = Human Muscle Cell Passage 3
  # 18 = Human Muscle Cell Passage 4
  # 19 = Human Mycoplasma
  # 20 = Human Muscle Histology
  # 21 = Human Adipose Histology
  
  # randomize batches for each assay & tissue
  # we want to randomize on site, randomgroupcode (which includes ped vs adult), sex, calculatedage 
  # all samples of a pid will stay together 
  # randomization will be independently performed in each tissue
  if(!"assay" %in% colnames(all_meta)){
    if(verbose){
      message("'assay' is not in the column names of the merged biospecimen and shipment metadata. Batching will assume that all samples from a given tissue are for a single assay.")
    }
    all_meta[,batching_group := sampletypecode]
  }else{
    all_meta[,batching_group := paste0(assay, '_', sampletypecode)]
  }
  if(length(unique(all_meta[,sampletypecode]))>1){
    if(verbose){
      message(sprintf("Batching will be performed separately for each of the following 'sampletypecode': %s\nIf this is not what you want, respecify 'api_metadata$SampleTypeCode' AND/OR 'shipment_manifest$`Sample Type`' columns in the input files.", paste0(unique(all_meta[,sampletypecode]), collapse=", ")))
    }
  }
  
  #### batching ##### 
  if(verbose) message()
  for (b in unique(all_meta[,batching_group])){
    
    if(verbose & length(unique(all_meta[,batching_group])) > 1){
      message(sprintf("\n\n--- BATCHING '%s' SAMPLES ---\n", b))
    }
    
    # check that the outfile hasn't already been generated
    outfile1 = sprintf("%s/files/%s-samples_UNBLINDED-batch-characteristics.csv", outdir, gsub(" ","-",b))
    if(file.exists(outfile1) & !overwrite){
      m = sprintf("File '%s' for '%s' samples already exists. Skipping. Set 'overwrite' to TRUE to ignore existing batching results and force rebalancing.",
                  outfile1, b)
      warning(m)
      next
    }
    
    # set up balancing strictness
    # can make this an argument in the future
    # reset this for each tissue 
    balance_strictness = rep(init_balance_strictness, length(balance_vars))
    names(balance_strictness) = balance_vars
    
    curr_batch = unique(all_meta[batching_group == b])
    curr_batch_pid = unique(curr_batch[,.(codedsiteid, pid, randomgroupcode, sex_psca, calculatedage, older_than_median)])
    
    # how many samples per person?
    curr_batch_n = data.table(table(curr_batch[,pid]))
    colnames(curr_batch_n) = c('pid','N')
    
    curr_batch_pid[,pid := as.character(pid)]
    curr_batch_n[,pid := as.character(pid)]
    curr_batch_pid = merge(curr_batch_pid, curr_batch_n, by='pid')
    
    # is it a single batch?
    if(sum(curr_batch_pid[,N])<=max_n_per_batch){
      message(sprintf("Number of '%s' samples (%s) is less than or equal to 'max_n_per_batch' (%s). No rebalancing necessary.",
                      b,
                      sum(curr_batch_pid[,N]),
                      max_n_per_batch))
      batches = curr_batch_pid
      batches[,batch := 1]
    }else{
      if(strict_size | max_full){
        batches = make_batches_strict(curr_batch_pid, max_n_per_batch, balance_strictness, max_full, 
                                      max_inner_loop_iter, max_outer_loop_iter, strict_size, 
                                      balance_vars = balance_vars,
                                      verbose = verbose)
      }else{
        batches = make_random_batches_not_strict(curr_batch_pid, max_n_per_batch, balance_strictness, b, 
                                                 max_full, max_inner_loop_iter, max_outer_loop_iter, 
                                                 balance_vars = balance_vars, verbose = verbose)
      }
    }
    
    # visualize batch balance
    if(all(c("calculatedage","sex_psca","codedsiteid","randomgroupcode") %in% colnames(batches))){
      # make a plot (specific to MoTrPAC)
      if(verbose){
        message("Inspect the heatmap saved to 'plots/' to assess batch balance. Rerun the function if you are not satisfied.")
      }
      p = plot_batches(batches, curr_batch, outdir, b)
      print(p)
    }else{
      if(verbose){
        message("Inspect the tables printed to the console to assess batch balance. Rerun the function if you are not satisfied.\nNote a heatmap is not generated because some expected columns are missing. Expected columns include:\n 'calculatedage','sex_psca','codedsiteid','randomgroupcode'.")
      }
      # print tables for specified balance vars 
      for(var in balance_vars){
        print(table(batches[,batch], batches[,get(var)], dnn = c("batch", var)))
      }
    }

    # write two versions to file 
    # batch characteristics 
    write.table(batches, file=outfile1, sep=',', col.names=T, row.names=F, quote=F)
    
    # current and new positions 
    curr_batch[,pid := as.character(pid)]
    all_info = merge(curr_batch, batches[,.(pid, batch)], by='pid')
    stopifnot(nrow(all_info) == nrow(curr_batch))
    
    if(block_randomization){
      all_info[,ptmp := NA_real_]
      for(BATCH in unique(all_info[,batch])){
        # within a batch, randomize patients while keeping all samples from a patient together
        pids = unique(all_info[batch==BATCH,pid])
        pids_rand = 1:length(pids)
        names(pids_rand) = as.character(sample(pids, length(pids), replace=F))
        all_info[batch==BATCH, ptmp := pids_rand[as.character(pid)]]
        # then randomize samples within patient
        rand_n = runif(nrow(all_info[batch==BATCH]))
        all_info[batch==BATCH, ptmp_rand := ptmp + rand_n]
      }
      all_info = all_info[order(batch, ptmp_rand, decreasing = F)]
      all_info[,injection_order := 1:nrow(all_info)]
      all_info[,c("ptmp","ptmp_rand") := NULL]
      positions = all_info[,.(injection_order, viallabel, barcode, sampletypecode, box, position, batch)]
      colnames(positions) = c('injection_order', 'viallabel','barcode','sampletypecode','shipping_box','shipping_position','new_batch')
    }else{
      positions = all_info[,.(viallabel, barcode, sampletypecode, box, position, batch)]
      colnames(positions) = c('viallabel','barcode','sampletypecode','shipping_box','shipping_position','new_batch')
    }
    positions[,new_batch := paste0("batch_",new_batch)]
    # order across rows 
    if(all(grepl("^[A-z]", positions[,shipping_position]))){
      positions[,shipping_row := sapply(shipping_position, function(x) unname(unlist(strsplit(x, '')))[1])] # first character
      positions[,shipping_column := sapply(shipping_position, function(x) as.numeric(paste(unname(unlist(strsplit(x, '')))[2:3],collapse='')))] # second and third characters
      positions = positions[order(shipping_box, shipping_row, shipping_column)]
      positions[,c('shipping_row','shipping_column') := NULL]
    }else if(all(grepl("^[0-9]", positions[,shipping_position]))){
      positions[,shipping_position := as.numeric(shipping_position)]
      positions = positions[order(shipping_box, shipping_position)]
    }else{
      warning("Shipping positions for '%s' samples are a combination of numbers and plate positions. How are we supposed to order these? %s",
              b,
              paste(unique(positions[,shipping_position], collapse=',')))
    }
    
    if(block_randomization){
      positions = positions[order(injection_order, decreasing=F)]
    }
    
    if(!separate_batch_files){
      write.table(positions, file=sprintf("%s/files/%s-samples_BLINDED-batch-assignments.csv", outdir, gsub(" ","-",b)), sep=',', col.names=T, row.names=F, quote=F)
    }else{
      for(batch in unique(positions[,new_batch])){
        sub = positions[new_batch == batch]
        write.table(sub, file=sprintf("%s/files/%s-samples_BLINDED-%s-assignments.csv", outdir, gsub(" ","-",b), batch), sep=',', col.names=T, row.names=F, quote=F)
      }
    }
  }
  message("Done!")
  message(sprintf("\nPlease manually check the plots in '%s/plots' to ensure that batches are satisfactorily balanced, i.e. numbers are reasonably distributed across each ROW/variable level. Rerun the function with 'overwrite=TRUE' if you are not satisfied with the balance.", gsub("/$","",outdir)))
}


#' Read in input files 
#'
#' Internal function used to read in .xlsx, .tsv, or .csv files. "Box"/"box" columns are coerced to character. 
#'
#' @param path character, path or URL to file, e.g. 'Stanford_ADU830-10060_120720.xlsx'
#'
#' @return data table
read_input = function(path){
  
  if(grepl("\\.xlsx$", path, ignore.case = TRUE)){
    data = data.table(read.xlsx(path, sheet = 1, sep.names = " ", check.names = FALSE))
  }else if(grepl("\\.tsv$", path, ignore.case = TRUE)){
    data = fread(file=path, sep="\t", header=TRUE)
  }else if(grepl("\\.csv$", path, ignore.case = TRUE)){
    data = fread(file=path, sep=",", header=TRUE)
  }else{
    stop(sprintf("'%s' is not a .xlsx, .tsv, or .csv file.", path))
  }
  if("box" %in% tolower(colnames(data))){
    whichcol = colnames(data)[which(tolower(colnames(data))=="box")]
    data[,(whichcol):= lapply(.SD, as.character), .SDcols = whichcol]
  }
  return(data)
}


#' Check batch balance
#' 
#' Internal function. 
#' Given a data table of batches, check if the batch is 
#' balanced given the \code{strictness} parameter. 
#'
#' @param curr_batch_pid data table of batches
#' @param strictness integer, initial strictness of balance checks, 
#'   with 10 being the strictest and 1 being the most lenient
#' @param max_full bool, whether to force as many batches as possible to have *exactly* \code{max_n_per_batch} samples
#' @param balance_vars character vector of variables for which to check balance 
#' @param verbose bool, whether to print progress messages
#'
#' @return a named list of two items:
#' \describe{
#'   \item{\code{success}}{bool, whether batches are balanced}
#'   \item{\code{failed}}{character vector of \code{balance_vars} that were not balanced at the specified strictness}
#' }
#'
#' @details TODO
check_batch_balance = function(curr_batch_pid, 
                               strictness, max_full, 
                               balance_vars = c('codedsiteid','randomgroupcode','sex_psca'),
                               verbose = TRUE){
  
  # for data.table variables 
  batch = NULL 
  
  leniency = lapply(strictness, function(x) max(10-x, 0))
  
  # is there only one batch?
  if(length(unique(curr_batch_pid[,batch]))==1){
    if(max_full){
      message("With 'max_full_batches' set to TRUE, batch balance is only checked in full batches because the remaining batch can have very few samples, in which case it would be difficult to assess balance. There is only one full batch for these samples, so batch balance will not be checked. Please inspect the generated plots to ensure that the batches are sufficiently balanced. If not, rerun the script until you are satisfied with the balance.")
      return(list(success = TRUE,
                  failed = c()))
    }else{
      warning("The code should never get here...")
    }
  }
  
  # lowest strictness 1 == just make sure more than one level per balance_vars are represented in each batch
  # this is the default for small batches
  # high strictness 10 == the optimal number of samples per level per balance_vars are represented in each batch
  # if necessary, iteratively lower strictness until a feasible scheme is found 
  
  redo = F
  failed = c()
  for(var in balance_vars){
    if(!var %in% colnames(curr_batch_pid)){
      warning(sprintf("Variable '%s' not found in column names of batch assignments. Skipping.", var))
      next
    }
    
    # balance
    m = as.matrix(table(curr_batch_pid[,get(var)], curr_batch_pid[,batch],  dnn=c("var","batch")))
    
    if(leniency[[var]] == 9){
      # just make sure at least two levels of each balance_var are present in each batch
      # does any batch have all the same value for one of these groups?
      m = as.matrix(table(curr_batch_pid[,get(var)], curr_batch_pid[,batch]))
      m[m > 0] = 1
      if(any(colSums(m) == 1)){
        redo = T
        failed = c(failed, var)
      }
    }else{
      # optimal balance per row/level
      optimal_per_row = apply(m, 1, function(x) floor(mean(x)))
      optimal_per_row = optimal_per_row - leniency[[var]]
      optimal_per_row[optimal_per_row < 0] = 0
      # check that every value per row has at least this value
      optimal_mat = matrix(data = optimal_per_row, nrow=nrow(m), ncol=ncol(m))
      # if there are more samples for a given level than 2xn_batches, require at least one sample per batch 
      optimal_mat[which(rowSums(m) > 2*ncol(m)),] = apply(optimal_mat[which(rowSums(m) > 2*ncol(m)),], c(1,2), function(x) max(x, 1))
      
      if(any(m < optimal_mat)){
        redo = T
        failed = c(failed, var)
      }
    }
  }
  
  if(redo){
    curr_batch_pid = NULL
  }else{
    if(verbose) message("Success!")
  }
  return(list(success = !redo,
              failed = failed))
}


#' ID optimal batch sizes
#'
#' Internal function. Given inputs, determine how many samples should be in each batch,
#' stipulating that all samples from the same patient stay in the same batch. 
#'
#' @param curr_batch_pid data table of input samples
#' @param max_n_per_batch integer, max number of samples per batch
#' @param max_full bool, whether to force as many batches as possible to have *exactly* \code{max_n_per_batch} samples
#'
#' @return list of optimal batch sizes
id_optimal_batch_sizes = function(curr_batch_pid, max_n_per_batch, max_full){
  
  # for data.table variables
  N = NULL
  
  optimal_n_batches = ceiling(sum(curr_batch_pid[,N]) / max_n_per_batch)
  if(!max_full){
    optimal_batch_sizes = list()
    for(i in 1:optimal_n_batches){
      optimal_batch_sizes[[i]] = 0
    }
    j = 1
    for(i in 1:sum(curr_batch_pid[,N])){
      optimal_batch_sizes[[j]] = optimal_batch_sizes[[j]] + 1
      j = j + 1
      if(j > optimal_n_batches){
        j = 1
      }
    }
    optimal_batch_sizes = unlist(optimal_batch_sizes)
  }else{
    # all but one batch should have max_n_per_batch
    n_full = floor(sum(curr_batch_pid[,N])/max_n_per_batch)
    remaining = sum(curr_batch_pid[,N]) - max_n_per_batch*n_full
    optimal_batch_sizes = c(rep(max_n_per_batch, n_full), remaining)
  }
  return(optimal_batch_sizes)
}


#' ID feasible batch sizes
#'
#' Internal function. 
#' If [id_optimal_batch_sizes()] fails, identify the minimum possible number of 
#' batches and their ideal sizes. Only used if \code{strict_size} is TRUE. 
#'
#' @param curr_batch_pid data table of input samples
#' @param max_n_per_batch integer, max number of samples per batch
#' @param max_full bool, whether to force as many batches as possible to have *exactly* \code{max_n_per_batch} samples
#'
#' @return vector of batch sizes
id_feasible_batch_sizes = function(curr_batch_pid, max_n_per_batch, max_full){
  
  # for data.table variables 
  N = pid = batch = NULL
  
  n_samples_per_pid = curr_batch_pid[,N]
  names(n_samples_per_pid) = curr_batch_pid[,pid]
  n_samples_per_pid = n_samples_per_pid[order(n_samples_per_pid, decreasing=T)]
  n_samples_per_pid = as.list(n_samples_per_pid)
  
  sets = list()
  i = 1
  while(length(n_samples_per_pid) > 0){
    new_set = n_samples_per_pid[1]
    n_samples_per_pid = n_samples_per_pid[-1] # remove from remaining samples
    while(sum(unlist(new_set)) != max_n_per_batch){
      needed = max_n_per_batch - sum(unlist(new_set))
      # find first sample with this many samples or less
      to_add = n_samples_per_pid[n_samples_per_pid <= needed][1]
      if(is.null(unlist(to_add))){
        sets[[i]] = new_set
        break
      }
      n_samples_per_pid[names(to_add)] = NULL
      new_set = c(new_set, to_add)
    }
    sets[[i]] = new_set
    i = i+1
  }
  
  # match pid to batch
  pid_to_batches = c()
  b2 = lapply(sets, unlist)
  for(i in 1:length(b2)){
    pids = names(b2[[i]])
    pid_to_batch = rep(i, length(pids))
    names(pid_to_batch) = pids
    pid_to_batches = c(pid_to_batches, pid_to_batch)
  }
  curr_batch_pid[,batch := pid_to_batches[pid]]
  
  sizes = unlist(lapply(sets, function(x) sum(unlist(x))))
  if(max_full){
    return(sizes)
  }
  # can we make the smallest batch less small?
  # what is the difference between the smallest batch and the next smallest batch?
  ordered_sizes = sizes[order(sizes, decreasing=F)]
  diff = ordered_sizes[2] - ordered_sizes[1]
  if(diff > 0){
    # take sample from a full batch and move to the small batch
    full_batches = which(sizes == max_n_per_batch)
    PID = curr_batch_pid[batch %in% full_batches & N <= diff, pid][1]
    pid_size = curr_batch_pid[pid==PID, N]
    # switch batch
    old_batch = curr_batch_pid[pid==PID, N]
    new_batch = which.min(sizes)
    curr_batch_pid[,batch := new_batch]
    # update batch size
    sizes[old_batch] = sizes[old_batch] - pid_size
    sizes[new_batch] = sizes[new_batch] + pid_size
  }
  
  # at this point we've IDed the min number of batches 
  # and feasible batch sizes
  return(sizes)
}


#' Make batches with strict size
#' 
#' Internal function. If possible, make all batches with exactly \code{max_n_per_batch} samples. 
#' Otherwise, opt for batches as close to the same size as possible. 
#'
#' @param curr_batch_pid data table of input samples
#' @param max_n_per_batch integer, max number of samples per batch
#' @param balance_strictness integer, initial strictness of balance checks, 
#'   with 10 being the strictest and 1 being the most lenient
#' @param max_full bool, whether to force as many batches as possible to have *exactly* \code{max_n_per_batch} samples
#' @param balance_vars character vector, force batches to include samples from at 
#'   least one level of each of these variables
#' @param verbose bool, whether to print progress messages
#' @param max_inner_loop_iter integer, max number of failed attempts to fit all samples 
#'   in batches before increasing the number of batches
#' @param max_outer_loop_iter integer, max number of failed attempts to find optimally 
#'   balanced batches before relaxing the stringency of the balance checks
#' @param strict_size bool, force *all* batches to be as close to \code{max_n_per_batch} as possible. 
#'   Most applicable for small batches (e.g. < 20).
#'
#' @return data table of batch assignments or NULL if batching failed 
#' 
#' @details TODO
make_batches_strict = function(curr_batch_pid, 
                               max_n_per_batch, 
                               balance_strictness, 
                               max_full,
                               max_inner_loop_iter, 
                               max_outer_loop_iter, 
                               strict_size, 
                               balance_vars = c('codedsiteid','randomgroupcode','sex_psca'),
                               verbose=TRUE){
  
  # for data.table variables 
  N = pid = batch = total = NULL
  
  feasible_batches = FALSE
  
  n_samples_per_pid = curr_batch_pid[,N]
  names(n_samples_per_pid) = curr_batch_pid[,pid]
  n_samples_per_pid = as.list(n_samples_per_pid)
  
  # first try to find batches with ideal batch size 
  batch_sizes = id_optimal_batch_sizes(curr_batch_pid, max_n_per_batch, max_full)
  
  # keep track of which variables are failing balancing
  if(max_full){
    failed_balancing = list()
    for(v in balance_vars){
      failed_balancing[[v]] = 0
    }
  }
  
  # is it possible to make batches of this size? give up after 1e6 tries
  outer_loop = 1
  inner_loop = 1
  balanced_batches = F
  while(!balanced_batches){
    
    # set up batches 
    new_sets = list()
    for(i in 1:length(batch_sizes)){
      new_sets[[i]] = list()
    }
    remaining_room_per_batch = batch_sizes
    names(remaining_room_per_batch) = 1:length(remaining_room_per_batch)
    
    # randomly reorder n_samples_per_pid
    n_samples_per_pid_reordered = n_samples_per_pid[sample(1:length(n_samples_per_pid), length(n_samples_per_pid), replace=F)]
    while(!all(remaining_room_per_batch==0)){ # while there are samples left to place...
      # extract the first set 
      set_to_place = n_samples_per_pid_reordered[1] 
      # put this in the first batch that it fits
      batch_to_place = as.numeric(names(remaining_room_per_batch)[remaining_room_per_batch>=set_to_place][1])
      if(is.na(batch_to_place)){
        break
      }
      new_sets[[batch_to_place]] = c(new_sets[[batch_to_place]], set_to_place)
      remaining_room_per_batch[[batch_to_place]] = remaining_room_per_batch[[batch_to_place]] - unlist(set_to_place)
      n_samples_per_pid_reordered = n_samples_per_pid_reordered[-1] # only remove the sample after it's been placed
    }
    # this worked? 
    if(all(remaining_room_per_batch==0)){
      # only print once in a while if inner_loop is small 
      if(inner_loop < 100){
        pmessage = outer_loop %% 100 == 0
      }else{
        pmessage = T
      }
      if(verbose & pmessage){
        message(sprintf("Identified the %sth combination of samples that fits the ideal batch sizes after %s iterations. Checking balance...", outer_loop, inner_loop))
      }
      # match pid to batch
      pid_to_batches = c()
      b2 = lapply(new_sets, unlist)
      for(i in 1:length(b2)){
        pids = names(b2[[i]])
        pid_to_batch = rep(i, length(pids))
        names(pid_to_batch) = pids
        pid_to_batches = c(pid_to_batches, pid_to_batch)
      }
      curr_batch_pid[,batch := pid_to_batches[pid]]
      
      # check if batches are balanced 
      if(max_full){
        # if max_full, ignore the smallest batch when checking batch balance
        curr_sizes = curr_batch_pid[,list(total=sum(N)), by=batch]
        smallest_batch = curr_sizes[which.min(total), batch]
        batches_to_check = curr_batch_pid[batch != smallest_batch]
      }else{
        batches_to_check = curr_batch_pid
      }
      check_batches = check_batch_balance(batches_to_check, balance_strictness, max_full, balance_vars, verbose)
      balanced_batches = check_batches$success
      batch_assignments = curr_batch_pid
      outer_loop = outer_loop + 1
      inner_loop = 1
      
      # keep track of which variables failed balancing (only for max_full)
      if(max_full){
        for(v in check_batches$failed){
          failed_balancing[[v]] = failed_balancing[[v]] + 1
        }
      }
    }else{
      inner_loop = inner_loop + 1
      if(inner_loop>max_inner_loop_iter & !feasible_batches){
        new_batch_sizes = id_feasible_batch_sizes(curr_batch_pid, max_n_per_batch, max_full)
        batch_sizes = new_batch_sizes
        warning(sprintf("With %s total samples and maximum %s samples per batch, the ideal number of samples per batch are as follows:\n%s\nHowever, after %s iterations, no combination of samples was found to fit these batch sizes. Trying again with the following batch sizes:\n%s",
                        sum(curr_batch_pid[,N]),
                        max_n_per_batch,
                        paste0(batch_sizes, collapse=', '),
                        inner_loop,
                        paste0(batch_sizes, collapse=', ')))
        feasible_batches = T
      }
    }
    
    # check if outer_loop is over the max
    if(outer_loop > max_outer_loop_iter){
      if(strict_size){
        message(sprintf("Target batch sizes: %s", paste0(batch_sizes, collapse=', ')))
        stop(sprintf("With %s total samples, maximum %s samples per batch, and target batch sizes printed above, well-balanced batches were not found in %s candidate batches. You are currently requiring all batches to have samples from more than one group for all of the following variables:\n    %s\nHope for better luck and rerun the script - OR - try removing the least important variable from this list using the 'vars_to_balance' argument and rerun the function.",
                     sum(curr_batch_pid[,N]),
                     max_n_per_batch,
                     outer_loop-1,
                     paste0(balance_vars, collapse=', ')))
      }
      if(max_full){
        too_strict = names(which.max(failed_balancing))
        message(sprintf("With %s total samples and %s samples per batch whenever possible, balanced batches were not found in %s candidate batches. The current balance strictness parameters are as follows:\n",
                        sum(curr_batch_pid[,N]),
                        max_n_per_batch,
                        outer_loop-1))
        print(balance_strictness)
        message(sprintf("Decreasing balance strictness for %s by 1.", too_strict))
        
        balance_strictness[[too_strict]] = balance_strictness[[too_strict]] - 1
        outer_loop = 1
        inner_loop = 1
        # reset table
        for(v in balance_vars){
          failed_balancing[[v]] = 0
        }
      }
    }
  }
  return(batch_assignments)
}


#' Make large batches with flexible size
#'
#' Internal function. TODO. 
#' 
#' @param curr_batch_pid data table of input samples
#' @param max_n_per_batch integer, max number of samples per batch
#' @param balance_strictness integer, initial strictness of balance checks, 
#'   with 10 being the strictest and 1 being the most lenient
#' @param b character, sample type description
#' @param balance_vars character vector, force batches to include samples from at 
#'   least one level of each of these variables
#' @param verbose bool, whether to print progress messages
#' @param max_full bool, whether to force as many batches as possible to have *exactly* \code{max_n_per_batch} samples
#' @param max_inner_loop_iter integer, max number of failed attempts to fit all samples 
#'   in batches before increasing the number of batches
#' @param max_outer_loop_iter integer, max number of failed attempts to find optimally 
#'   balanced batches before relaxing the stringency of the balance checks
#'
#' @return data table of batch assignments or NULL if batching failed 
#'
#' @details TODO
make_random_batches_not_strict = function(curr_batch_pid, max_n_per_batch, balance_strictness, b, 
                                          max_full, max_inner_loop_iter, max_outer_loop_iter, 
                                          balance_vars = c('codedsiteid','randomgroupcode','sex_psca'),
                                          verbose=T){
  
  # data.table variables 
  N = pid = batch = NULL
  
  n_samples_per_pid = curr_batch_pid[,N]
  names(n_samples_per_pid) = curr_batch_pid[,pid]
  n_samples_per_pid = as.list(n_samples_per_pid)
  n_batches = ceiling(sum(curr_batch_pid[,N]) / max_n_per_batch)
  
  # keep track of which variables are failing balancing
  failed_balancing = list()
  for(v in balance_vars){
    failed_balancing[[v]] = 0
  }
  
  outer_loop = 1
  inner_loop = 1
  balanced_batches = F
  already_reduced_stringency = F
  while(!balanced_batches){
    
    # set up batches 
    new_sets = list()
    for(i in 1:n_batches){
      new_sets[[i]] = list()
    }
    remaining_room_per_batch = rep(max_n_per_batch, n_batches)
    names(remaining_room_per_batch) = 1:length(remaining_room_per_batch)
    
    # randomly reorder n_samples_per_pid
    n_samples_per_pid_reordered = n_samples_per_pid[sample(1:length(n_samples_per_pid), length(n_samples_per_pid), replace=F)]
    batch_iter = 1
    while(length(n_samples_per_pid_reordered)>0){ # while there are samples left to place...
      # extract the first set 
      set_to_place = n_samples_per_pid_reordered[1] 
      # put this in the first batch that it fits
      which_fits = unname(which(remaining_room_per_batch>=set_to_place))
      if(length(which_fits)==0){
        break
      }
      if(any(which_fits >= batch_iter)){
        # get first one greater than or equal to batch_iter
        which_fits = which_fits[which_fits >= batch_iter]
        which_fits = min(which_fits)
      }else{ 
        # if that doesn't work, start again at 1 
        which_fits[which_fits < batch_iter]
        which_fits = min(which_fits)
      }
      new_sets[[which_fits]] = c(new_sets[[which_fits]], set_to_place)
      remaining_room_per_batch[[which_fits]] = remaining_room_per_batch[[which_fits]] - unlist(set_to_place)
      n_samples_per_pid_reordered = n_samples_per_pid_reordered[-1]  # only remove the sample after it's been placed
      batch_iter = batch_iter + 1
      if(batch_iter > n_batches){
        batch_iter = 1
      }
    }
    # this worked? 
    if(length(n_samples_per_pid_reordered)==0){
      # only print once in a while if inner_loop is small 
      if(inner_loop < 100){
        pmessage = outer_loop %% 500 == 0
      }else{
        pmessage = T
      }
      if(verbose & pmessage){
        message(sprintf("Identified the %sth combination of samples that fits the ideal batch sizes after %s iterations. Checking balance...", outer_loop, inner_loop))
      }
      
      if(inner_loop > 5000 & !already_reduced_stringency){
        message(sprintf("It took at least 5000 iterations to find a *single* combination of '%s' samples that fits in the ideal number of batches. Reducing stringency for batch balance checks.", b))
        for(v in names(balance_strictness)){
          balance_strictness[[v]] = 1
        }
        already_reduced_stringency = T
      }
      
      # match pid to batch
      pid_to_batches = c()
      b2 = lapply(new_sets, unlist)
      for(i in 1:length(b2)){
        pids = names(b2[[i]])
        pid_to_batch = rep(i, length(pids))
        names(pid_to_batch) = pids
        pid_to_batches = c(pid_to_batches, pid_to_batch)
      }
      stopifnot(length(pid_to_batches) == nrow(curr_batch_pid))
      curr_batch_pid[,batch := pid_to_batches[pid]]
      
      # check if batches are balanced 
      check_batches = check_batch_balance(curr_batch_pid, balance_strictness, max_full, balance_vars, verbose)
      balanced_batches = check_batches$success
      #batch_assignments = check_batches$batch_assignments
      batch_assignments = curr_batch_pid
      outer_loop = outer_loop + 1
      inner_loop = 1
      
      # keep track of which variables failed balancing 
      for(v in check_batches$failed){
        failed_balancing[[v]] = failed_balancing[[v]] + 1
      }
    }else{
      inner_loop = inner_loop + 1
      if(inner_loop > max_inner_loop_iter){
        warning("Can't find a combination of samples that fits in this many batches. Increasing the number of batches by 1.")
        n_batches = n_batches + 1
      }
    }
    
    # adjust balance strictness if necessary, and reset outer loop counter 
    if(outer_loop > max_outer_loop_iter){
      if(max_n_per_batch < 20){
        stop(sprintf("It looks like you want small batch sizes (N = %s). Please try re-running the script with strict_size = TRUE for a batching method more suitable for small batch sizes.",
                     max_n_per_batch))
      }else{
        too_strict = names(which.max(failed_balancing))
        message(sprintf("With %s total samples and up to %s samples per batch, balanced batches were not found in %s candidate batches. The current balance strictness parameters are as follows:\n",
                        sum(curr_batch_pid[,N]),
                        max_n_per_batch,
                        outer_loop-1))
        print(balance_strictness)
        if(all(balance_strictness == 10)){
          message(sprintf("Decreasing balance strictness for %s by 1.", paste0(names(balance_strictness), collapse=', ')))
          for(v in names(balance_strictness)){
            balance_strictness[[v]] = max(1, balance_strictness[[v]] - 1)
          }
        }else{
          message(sprintf("Decreasing balance strictness for %s by 1.", too_strict))
          balance_strictness[[too_strict]] = max(1, balance_strictness[[too_strict]] - 1)
        }
        
        outer_loop = 1
        inner_loop = 1
        # reset table
        for(v in balance_vars){
          failed_balancing[[v]] = 0
        }
      }
    }
    
  }
  return(batch_assignments)
}


#' Plot batch characteristics
#' 
#' Plot a heatmap of the numbers of samples per level per balance variable across batches. 
#' Well-balanced batches are indicated by relatively similar numbers across each ROW of the plot. 
#' This function saves the plot locally and returns the [pheatmap::pheatmap()] object. 
#'
#' @param batches batching output
#' @param curr_batch batching input
#' @param outdir character, path to output directory
#' @param b character, name of samples being batched
#' 
#' @seealso [batch_samples()]
#'
#' @return pheatmap object
#' 
#' @details 
#' Required columns in \code{batches} include 'calculatedage','sex_psca','codedsiteid','randomgroupcode','batch','N'. 
#' Required columns in \code{curr_batch} include 'codedsiteid','randomgroupcode'. 
plot_batches = function(batches, curr_batch, outdir, b){
  
  # make nice heatmap table
  b2 = copy(batches)
  # order batch numerically
  nb = length(unique(b2[,batch]))
  b2[, batch := factor(paste0("Batch ", batch), levels = c(paste0("Batch ", 1:nb)))]
  b2[,subj := 1]
  b2_tibble = gtsummary::as_tibble(b2[,.(calculatedage, sex_psca, codedsiteid, randomgroupcode, batch, N, subj)])
  tb = gtsummary::tbl_summary(b2_tibble, 
                              by = 'batch',
                              type = list(N ~ "continuous",
                                          calculatedage ~ "continuous",
                                          sex_psca ~ "categorical",
                                          codedsiteid ~ "categorical",
                                          randomgroupcode ~ "categorical",
                                          subj ~ "continuous"),
                              label = list(N ~ "N samples",
                                           calculatedage ~ "Age",
                                           sex_psca ~ "Sex",
                                           codedsiteid ~ "Site code",
                                           randomgroupcode ~ "Intervention group",
                                           subj ~ "N subjects"),
                              statistic = list(all_continuous() ~ "{median} ({min},{max})",
                                               all_categorical() ~ "{n}",
                                               N ~ "{sum}",
                                               subj ~ "{sum}"),
                              digits = gtsummary::all_continuous() ~ 0)
  tb_df = as.data.frame(gtsummary::as_tibble(tb))
  colnames(tb_df) = c('characteristic', paste0('Batch', 1:nb))
  # can't figure out how to color this or save it to a file...
  # can we just use pheatmap for now?
  rownames(tb_df) = tb_df$characteristic
  tb_df$characteristic = NULL
  
  # reorder rows
  sites = unique(curr_batch[,codedsiteid])
  sites = as.character(sites[order(sites, decreasing=F)])
  randgroup = unique(curr_batch[,randomgroupcode])
  randgroup = randgroup[order(randgroup)]
  tb_df = tb_df[c("N subjects", "N samples", "Age", "Sex", "1", "2", 
                  "Site code", sites,
                  "Intervention group", randgroup),]
  
  labels = tb_df
  labels[is.na(labels)] = ''
  values = tb_df
  values["Age",] = as.numeric(gsub(" .*","",values["Age",]))
  values = as.data.frame(apply(values, c(1,2), as.numeric))
  
  rownames(values)[rownames(values)=='Age'] = 'Age [med (min,max)]'
  rownames(labels)[rownames(labels)=='Age'] = 'Age [med (min,max)]'
  
  p = pheatmap(values,
               color = rev(heat.colors(50)),
               scale = "row",
               cluster_rows = F,
               cluster_cols = F,
               legend = F,
               display_numbers = labels,
               na_col = 'gray',
               fontsize_col = 10,
               fontsize_row = 10,
               fontsize_number = 10,
               main = b,
               gaps_row = 2)
  
  w = 0.5 + length(unique(batches[,batch]))*0.9
  pdf(sprintf("%s/plots/batch-characteristics_%s.pdf", outdir, b), width=w, height=8)
  print(p)
  dev.off()
  
  return(p)
}
