#!/bin/R
# Nicole Gay
# 22 April 2021
# updated 4 August 2021
# batching for MoTrPAC clinical samples

library(data.table)
library(readxl)
library(testit)
library(argparse)
library(ggplot2)
library(gtsummary)
library(pheatmap)

# ARGUMENTS #########################################################################################################

### if you want to run this code in RStudio instead of from the command line, 
### comment out this chunk and uncomment/define the variables listed underneath
parser = ArgumentParser()
parser$add_argument("-ship", "--shipment-manifest-excel", required=T, type="character", nargs="+",
                    help="Path(s) to shipment manifest Excel files, e.g. Stanford_ADU830-10060_120720.xlsx Stanford_PED830-10062_120720.xlsx")
parser$add_argument("-api", "--api-metadata-csv", required=T, type="character", nargs="+",
                    help="Path(s) to sample metadata from web API, e.g. ADU830-10060.csv PED830-10062.csv")
parser$add_argument("-max", "--max-n-per-batch", type="integer", required=T,
                    help="Max number of samples per batch")
parser$add_argument("-s", "--strict-size", action="store_true", default=F,
                    help="Force *all* batches to be as close to --max-n-per-batch as possible. Most applicable for small batches (e.g. < 20)")
parser$add_argument("-f", "--max-full-batches", action="store_true", default=F,
                    help="Force as many batches as possible to have *exactly* --max-n-per-batch samples.")
parser$add_argument("-v", "--vars-to-balance", type="character", default=c('codedsiteid','randomgroupcode','sex_psca','older_than_40'),
                    help="Force batches to include samples from at least two groups of each of these variables. Must be defined in --api-metadata-csv")
parser$add_argument("-o", "--outdir", type="character", default=".",
                    help="Path to output directory")
parser$add_argument("-q", "--quietly", action="store_true", default=FALSE,
                    help="Silence progress messages")
parser$add_argument("-inner", "--max-inner-loop-iter", type="integer", default=1e6,
                    help="Max number of failed attempts to fit all samples in batches before increasing the number of batches")
parser$add_argument("-outer", "--max-outer-loop-iter", type="integer", default=NULL,
                    help="Max number of failed attempts to find optimally balanced bacthes before relaxing the stringency of the balance checks")
parser$add_argument("-b", "--balance-strictness", type="integer", default=NULL,
                    help="Initial strictness of balance checks, with 10 being the strictest and 1 being the most lenient")
parser$add_argument("--overwrite", action="store_true", default=FALSE,
                    help="Overwrite existing batching results")
parser$add_argument("--tissue-subset", default=NULL,
                    help="Run batching for a single tissue. Must be a value in the 'Sample Type' column of one --shipment-manifest-excel OR a value in the 'SampleTypeCode' column of one --api-metadata-csv")
parser$add_argument("--block-randomization", action="store_true", default=FALSE,
                    help="Block randomization for metabolomics sites: samples within a batch are ordered by individual; samples within an individual are randomized. This adds an 'injection_order' column.")
parser$add_argument("--separate-batch-files", action="store_true", default=FALSE,
                    help="Write separate BLINDED output files per batch")

args = parser$parse_args()
shipments = args$shipment_manifest_excel
apis = args$api_metadata_csv
max_n_per_batch = args$max_n_per_batch
strict_size = args$strict_size
max_full = args$max_full_batches
balance_vars = args$vars_to_balance
outdir = args$outdir
verbose = !args$quietly
max_inner_loop_iter = args$max_inner_loop_iter
max_outer_loop_iter = args$max_outer_loop_iter
overwrite = args$overwrite
init_balance_strictness = args$balance_strictness
tissue_subset = args$tissue_subset
block_randomization = args$block_randomization
separate_batch_files = args$separate_batch_files

####
#
# # Broad proteomics example
# shipments = "~/Desktop/broad_batches/ShipmentContents_BroadCarr_012521.xlsx"
# apis = "~/Desktop/broad_batches/ADU822-10074.csv"
# max_n_per_batch = 15
# strict_size = T
# max_full = F
# balance_vars = c('codedsiteid','randomgroupcode','sex_psca','older_than_40')
# outdir = "~/Desktop/broad_batches"
# verbose = T
# max_inner_loop_iter = 1e6
# max_outer_loop_iter = 1000
# overwrite = F
# init_balance_strictness = 1
# tissue_subset = NULL
# block_randomization = F
# separate_batch_files = F
#
# # Stanford GET example
# shipments = c("/Users/nicolegay/Documents/motrpac/CLINICAL/metadata/shipment/Stanford_ADU830-10060_120720.xlsx",
#               "/Users/nicolegay/Documents/motrpac/CLINICAL/metadata/shipment/Stanford_PED830-10062_120720.xlsx")
# apis = c("/Users/nicolegay/Documents/motrpac/CLINICAL/metadata/manifests/ADU830-10060.csv",
#          "/Users/nicolegay/Documents/motrpac/CLINICAL/metadata/manifests/PED830-10062.csv")
# max_n_per_batch = 94
# strict_size = F
# max_full = F
# balance_vars = c('codedsiteid','randomgroupcode','sex_psca','older_than_40')
# outdir = "~/Desktop/stanford_batches"
# verbose = T
# max_inner_loop_iter = 1e6
# max_outer_loop_iter = 5000
# overwrite = F
# init_balance_strictness = 8
# tissue_subset = 'Muscle'
# block_randomization = T
# separate_batch_files = T
# 
# # metab example
# shipments = c("/Users/nicolegay/Documents/motrpac/CLINICAL/metadata/shipment/Stanford_ADU830-10060_120720.xlsx")
# apis = c("/Users/nicolegay/Documents/motrpac/CLINICAL/metadata/manifests/ADU830-10060.csv")
# max_n_per_batch = 88
# strict_size = F
# max_full = T
# balance_vars = c('codedsiteid','randomgroupcode','sex_psca','older_than_40')
# outdir = "~/Desktop/metab_batches"
# verbose = T
# max_inner_loop_iter = 1e6
# max_outer_loop_iter = 5000
# overwrite = F
# init_balance_strictness = 10
# tissue_subset = NULL
# block_randomization = F
# separate_batch_files = F
#
#### 

if(is.null(init_balance_strictness)){
  if(strict_size){
    init_balance_strictness = 1
  }else{
    init_balance_strictness = 10
  }
}

if(max_full & strict_size){
  stop("This script is not currently designed for cases where *both* --max-full-batches and --strict-size are TRUE. Please use only one of the two flags. --strict-size should be used for small batches, e.g. --max-n-per-batch < 20; --max-full-batches should be used for larger batches where the user wants to maximize the number of batches with exactly --max-n-per-batch samples.")
}

if(is.null(max_outer_loop_iter)){
  if(strict_size){
    max_outer_loop_iter = 1000
  }else{
    max_outer_loop_iter = 5000 
  }
}

# FUNCTIONS #########################################################################################################

# read in shipment manifests
# error reading in Excel file because "Box" is a mix of numeric and string. fix that
read_shipment = function(path){
  # read first line
  trunc = read_excel(path, n_max=2)
  which(colnames(trunc)=='Box')
  types = rep("guess", ncol(trunc))
  types[which(colnames(trunc)=='Box')] = 'text'
  # make sure "Box" is read in as character
  full = data.table(read_excel(path, col_types = types))
  return(full)
}


check_batch_balance = function(curr_batch_pid, 
                               strictness, 
                               balance_vars = c('codedsiteid','randomgroupcode','sex_psca')){
  
  leniency = lapply(strictness, function(x) max(10-x, 0))
  
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
    m = as.matrix(table(curr_batch_pid[,get(var)], curr_batch_pid[,batch]))
    
    if(leniency[[var]] == 9){
      # just make sure at least two levels of each balance_var are present in each batch
      # does any batch have all the same value for one of these groups?
      m = as.matrix(table(curr_batch_pid[,get(var)], curr_batch_pid[,batch]))
      m[m > 0] = 1
      if(any(colSums(m) == 1)){
        # if(verbose){
        #   message(sprintf("Batch(es) %s only includes samples from one group of %s:",
        #                   paste0(which(colSums(m) == 1), collapse=','),
        #                   var))
        #   print(table(curr_batch_pid[,get(var)], curr_batch_pid[,batch]))
        # }
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
        # if(verbose){
        #   message(sprintf("Batches not optimally balanced for %s with a strictness of %s", var, strictness))
        # }
        redo = T
        failed = c(failed, var)
      }
    }
  }
  
  if(redo){
    #if(verbose) message("Trying again to find more balanced batches...\n")
    curr_batch_pid = NULL
  }else{
    if(verbose) message("Success!")
  }
  return(list(success = !redo,
              failed = failed))
}


id_optimal_batch_sizes = function(curr_batch_pid, max_n_per_batch, max_full){
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


id_feasible_batch_sizes = function(curr_batch_pid, b, max_n_per_batch, max_full){
  
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


# minimize the number of batches
# when N samples is not divisible by max_n_per_batch, 
# opt for batches as close to the same size as possible 
make_batches_strict = function(curr_batch_pid, b, max_n_per_batch, balance_strictness, max_full,
                               balance_vars = c('codedsiteid','randomgroupcode','sex_psca'),
                               verbose=T){
  
  feasible_batches = F
  
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
      check_batches = check_batch_balance(batches_to_check, balance_strictness, balance_vars)
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
        new_batch_sizes = id_feasible_batch_sizes(curr_batch_pid, b, max_n_per_batch, max_full)
        batch_sizes = new_batch_sizes
        warning(sprintf("With %s total samples and maximum %s samples per batch, the ideal number of samples per batch are as follows:\n%s\nHowever, after %s iterations, no combination of samples was found to fit these batch sizes. Trying again with the following batch sizes:\n%s",
                        sum(curr_batch_pid[,N]),
                        max_n_per_batch,
                        paste0(batch_sizes, collapse=', '),
                        inner_loop,
                        paste0(batch_sizes, collapse=', ')))
        feasible_batches = T
      }
      if(outer_loop > max_outer_loop_iter){
        if(verbose){
          writeLines(paste0(batch_sizes, collapse=', '))
        }
        
        if(strict_size){
          stop(sprintf("With %s total samples, maximum %s samples per batch, and target batch sizes printed above, well-balanced batches were not found in %s candidate batches. You are currently requiring all batches to have samples from more than one group for all of the following variables:\n    %s\nHope for better luck and rerun the script - OR - try removing the least important variable from this list using the --vars-to-balance flag and rerun the script.",
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
  }
  return(batch_assignments)
}


make_random_batches_not_strict = function(curr_batch_pid, b, max_n_per_batch, balance_strictness,
                                          balance_vars = c('codedsiteid','randomgroupcode','sex_psca'),
                                          verbose=T){
  
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
        pmessage = outer_loop %% 100 == 0
      }else{
        pmessage = T
      }
      if(verbose & pmessage){
        message(sprintf("Identified the %sth combination of samples that fits the ideal batch sizes after %s iterations. Checking balance...", outer_loop, inner_loop))
      }
      
      if(inner_loop > 5000 & !already_reduced_stringency){
        message(sprintf("It took at least 1000 iterations to find a *single* combination of '%s' samples that fits in the ideal number of batches. Reducing stringency for batch balance checks.", b))
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
      assert(length(pid_to_batches) == nrow(curr_batch_pid))
      curr_batch_pid[,batch := pid_to_batches[pid]]
      
      # check if batches are balanced 
      check_batches = check_batch_balance(curr_batch_pid, balance_strictness, balance_vars)
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
      if(outer_loop > max_outer_loop_iter){
        if(max_n_per_batch < 20){
          stop(sprintf("It looks like you want small batch sizes (N = %s). Please try re-running the script with the --strict-size flag for a batching method more suitable for small batch sizes.",
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
  }
  return(batch_assignments)
}

# CHECK FORMATS #########################################################################################################

if(!all(sapply(shipments, function(x) grepl("\\.xls", x, ignore.case=T)))){
  stop(sprintf("Shipment manifests are not in the expected .xls or .xlsx format: %s", paste(shipments, collapse=', ')))
}
if(!all(sapply(apis, function(x) grepl("\\.csv$", x, ignore.case = T)))){
  stop(sprintf("API metadata files are not in the expected .csv format: %s", paste(apis, collapse=', ')))
}

# LOAD MANIFESTS #########################################################################################################

ship_list = list()
for (s in shipments){
  ship_list[[s]] = read_shipment(s)
}
ship = rbindlist(ship_list, fill=T)

# read in API metadata 
api_list = list()
for (a in apis){
  api_list[[a]] = fread(a, sep=',', header=T)
}
api = rbindlist(api_list, fill=T)

# make colnames lowercase for simplicity
colnames(ship) = tolower(colnames(ship))
colnames(api) = tolower(colnames(api))

api[,older_than_40 := calculatedage > 40]

# merge
api[,viallabel := as.character(viallabel)]
ship[,viallabel := as.character(viallabel)]
ship[,`2d barcode` := as.character(`2d barcode`)]
# PBMCs don't have barcodes?
all_meta = merge(api, ship, by='viallabel')
assert(nrow(all_meta) == nrow(ship))
all_meta[all_meta=='.'] = NA
               
if(nrow(all_meta)!=nrow(api)){
  warning(sprintf("The number of samples in the merged metadata (%s) is not the same as the number of samples in the biospecimen metadata '%s' (%s). Check for a merging error.\n",
          nrow(all_meta),
          paste0(basename(apis), collapse=', '),
          nrow(api)))
}

# check if we need to subset by a tissue
if(!is.null(tissue_subset)){
  # try to find a matching value in either `sample type` or `sampletypecode`
  if(tissue_subset %in% all_meta[,sampletypecode]){
    all_meta = all_meta[sampletypecode == tissue_subset]
  }else if(tissue_subset %in% all_meta[,`sample type`]){
    all_meta = all_meta[`sample type` == tissue_subset]
  }else if(tolower(tissue_subset) %in% tolower(all_meta[,sampletypecode])){
    all_meta = all_meta[tolower(sampletypecode) == tolower(tissue_subset)]
  }else if(tolower(tissue_subset) %in% tolower(all_meta[,`sample type`])){
    all_meta = all_meta[tolower(`sample type`) == tolower(tissue_subset)]
  }else if(is.numeric(all_meta[,`sample type`])){
    if(as.numeric(tissue_subset) %in% all_meta[,`sample type`]){
      all_meta = all_meta[`sample type` == as.numeric(tissue_subset)]
    }
  }else if(is.numeric(all_meta[,sampletypecode])){
    if(as.numeric(tissue_subset) %in% all_meta[,sampletypecode]){
      all_meta = all_meta[sampletypecode == as.numeric(tissue_subset)]
    }
  }else{
    stop(sprintf("Unable to match --tissue-subset '%s' to a value in either api$SampleTypeCode or ship$`Sample Type`. Available options to subset by tissue: %s",
                 tissue_subset, 
                 paste0(c(unique(all_meta[,sampletypecode]), unique(all_meta[,`sample type`])), collapse=', ')))
  }
}

# subset to existing samples
all_meta = all_meta[!(is.na(viallabel) | is.na(box))]

# SETUP #########################################################################################################

# make outdirs
system(sprintf("mkdir -p %s/plots", outdir))
system(sprintf("mkdir -p %s/files", outdir))

# "study" is redundant with "randomgroupcode"
# there can be multiple bid per pid. pid id human participant id. use pid as identifier
# visitcode = baseline versus post. ignore this because all samples from an individual will be together

if("assay" %in% colnames(all_meta)){
  all_meta = all_meta[,.(viallabel, pid, protocol, codedsiteid, barcode, 
                         sampletypecode, randomgroupcode, sex_psca, calculatedage, older_than_40, 
                         box, position, assay)]
}else{
  all_meta = all_meta[,.(viallabel, pid, protocol, codedsiteid, barcode, 
                         sampletypecode, randomgroupcode, sex_psca, calculatedage, older_than_40, 
                         box, position)]
}

# columns with 0 variance
remove = c()
for(c in colnames(all_meta)){
  if(length(unique(all_meta[,get(c)]))==1){
    remove = c(remove, c)
  }
}
message(sprintf("\nThese columns have 0 variance in the merged metadata: %s",
                paste0(remove, collapse=', ')))
# all_meta[,(remove) := NULL]

#table(all_meta[,sampletypecode], all_meta[,randomgroupcode])
# 4 = PaxGene RNA
# 5 = PBMC
# 6 = Muscle

# randomize batches for each assay & tissue
# we want to randomize on site, randomgroupcode (which includes ped vs adult), sex, calculatedage 
# all samples of a pid will stay together 
# randomization will be independently performed in each tissue
if(!"assay" %in% colnames(all_meta)){
  message("\n'assay' is not in the column names of the merged biospecimen and shipment metadata. Batching will assume that all samples from a given tissue are for a single assay.")
  all_meta[,batching_group := sampletypecode]
}else{
  all_meta[,batching_group := paste0(assay, '_', sampletypecode)]
}

# BATCHING #########################################################################################################

for (b in unique(all_meta[,batching_group])){

  if(verbose){
    message(sprintf("\n\n--- BATCHING '%s' SAMPLES ---\n", b))
  }

  # check that the outfile hasn't already been generated
  outfile1 = sprintf("%s/files/precovid_%s-samples_UNBLINDED-batch-characteristics.csv", outdir, gsub(" ","-",b))
  if(file.exists(outfile1) & !overwrite){
    m = sprintf("File %s for '%s' samples already exists. Skipping. Use the --overwrite flag to ignore existing batching results and force rebalancing.",
                outfile1, b)
    message(m)
    warning(m)
    next
  }

  # set up balancing strictness
  # can make this an argument in the future
  # reset this for each tissue 
  balance_strictness = rep(init_balance_strictness, length(balance_vars))
  names(balance_strictness) = balance_vars

  curr_batch = unique(all_meta[batching_group == b])
  curr_batch_pid = unique(curr_batch[,.(codedsiteid, pid, randomgroupcode, sex_psca, calculatedage, older_than_40)])
  
  # how many samples per person?
  curr_batch_n = data.table(table(curr_batch[,pid]))
  colnames(curr_batch_n) = c('pid','N')

  curr_batch_pid[,pid := as.character(pid)]
  curr_batch_n[,pid := as.character(pid)]
  curr_batch_pid = merge(curr_batch_pid, curr_batch_n, by='pid')
  
  if(strict_size | max_full){
    batches = make_batches_strict(curr_batch_pid, b, max_n_per_batch, balance_strictness, max_full, 
                                  balance_vars = balance_vars,
                                  verbose = verbose)
  }else{
    batches = make_random_batches_not_strict(curr_batch_pid, b, max_n_per_batch, balance_strictness,
                                             balance_vars = balance_vars, verbose = verbose)
  }

  # make nice heatmap table
  b2 = copy(batches)
  b2[, batch := paste0("Batch ", batch)]
  b2[,subj := 1]
  tb = gtsummary::tbl_summary(b2[,.(calculatedage, sex_psca, codedsiteid, randomgroupcode, batch, N, subj)], 
                         by='batch',
                         type=list(N ~ "continuous",
                                   calculatedage ~ "continuous",
                                   sex_psca ~ "categorical",
                                   codedsiteid ~ "categorical",
                                   randomgroupcode ~ "categorical",
                                   subj ~ "continuous"),
                         label=list(N ~ "N samples",
                                   calculatedage ~ "Age",
                                   sex_psca ~ "Sex",
                                   codedsiteid ~ "Site code",
                                   randomgroupcode ~ "Intervention group",
                                   subj ~ "N subjects"),
                         statistic = list(all_continuous() ~ "{median} ({min},{max})",
                                          all_categorical() ~ "{n}",
                                          N ~ "{sum}",
                                          subj ~ "{sum}"),
                         digits = list(subj ~ "0"))
  tb_df = as.data.frame(gtsummary::as_tibble(tb))
  colnames(tb_df) = c('characteristic', paste0('Batch', 1:max(batches[,batch])))
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
  
  # write two versions to file 
  # batch characteristics 
  write.table(batches, file=outfile1, sep=',', col.names=T, row.names=F, quote=F)
  
  # current and new positions 
  curr_batch[,pid := as.character(pid)]
  all_info = merge(curr_batch, batches[,.(pid, batch)], by='pid')
  assert(nrow(all_info) == nrow(curr_batch))
  
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
    write.table(positions, file=sprintf("%s/files/precovid_%s-samples_BLINDED-batch-assignments.csv", outdir, gsub(" ","-",b)), sep=',', col.names=T, row.names=F, quote=F)
  }else{
    for(batch in unique(positions[,new_batch])){
      sub = positions[new_batch == batch]
      write.table(sub, file=sprintf("%s/files/precovid_%s-samples_BLINDED-%s-assignments.csv", outdir, gsub(" ","-",b), batch), sep=',', col.names=T, row.names=F, quote=F)
    }
  }
}
message("\nDone!")
message(sprintf("\nPlease manually check the plots in %s/plots to ensure that batches are satisfactorily balanced, i.e. numbers are reasonably distributed across each ROW/variable level. Rerun the script if you are not satisfied with the balance.\n", gsub("/$","",outdir)))

warnings()
