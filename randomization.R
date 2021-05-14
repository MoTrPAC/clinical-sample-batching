#!/bin/R
# Nicole Gay
# 22 April 2021
# batching for MoTrPAC clinical samples

library(data.table)
library(readxl)
library(testit)
library(argparse)

# ARGUMENTS #########################################################################################################

### if you want to run this code in RStudio instead of from the command line, 
### comment out this chunk and uncomment/define the variables listed underneath
# parser = ArgumentParser()
# parser$add_argument("-s", "--shipment-manifest-excel", required=T, type="character", nargs="+",
#                     help="Path(s) to shipment manifest Excel files, e.g. Stanford_ADU830-10060_120720.xlsx Stanford_PED830-10062_120720.xlsx")
# parser$add_argument("-a", "--api-metadata-csv", required=T, type="character", nargs="+",
#                     help="Path(s) to sample metadata from web API, e.g. ADU830-10060.csv PED830-10062.csv")
# parser$add_argument("-n", "--max-n-per-batch", type="integer", required=T,
#                     help="Max number of samples per batch")
# parser$add_argument("-t", "--strict-size", action="store_true", default=F,
#                     help="Force all batches to be as close to --max-n-per-batch as possible. Most applicable for small batches (e.g. < 20)")
# parser$add_argument("-b", "--vars-to-balance", type="character", default=c('codedsiteid','randomgroupcode','sex_psca','older_than_40'),
#                     help="Force batches to include samples from at least two groups of each of these variables. Must be defined in --api-metadata-csv")
# parser$add_argument("-o", "--outdir", type="character", default=".",
#                     help="Path to output directory")
# parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
#                     help="Print progress messages [default]")
# args = parser$parse_args()
# shipments = args$shipment_manifest_excel
# apis = args$api_metadata_csv
# max_n_per_batch = args$max_n_per_batch
# strict_size = args$strict_size 
# balance_vars = args$vars_to_balance
# outdir = args$outdir
# verbose = args$verbose
# 
# 
###
shipments = "~/Desktop/broad_batches/ShipmentContents_BroadCarr_012521.xlsx"
apis = "~/Desktop/broad_batches/ADU822-10074.csv"
max_n_per_batch = 15
strict_size = T
balance_vars = c('codedsiteid','randomgroupcode','sex_psca','older_than_40')
outdir = "~/Desktop/broad_batches"
verbose = F
### 

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
                               balance_vars = c('codedsiteid','randomgroupcode','sex_psca')){
  redo = F
  for(var in balance_vars){
    if(!var %in% colnames(curr_batch_pid)){
      warning(sprintf("Variable '%s' not found in column names of batch assignments. Skipping.", var))
      next
    }
    # does any batch have all the same value for one of these groups?
    m = as.matrix(table(curr_batch_pid[,get(var)], curr_batch_pid[,batch]))
    m[m > 0] = 1
    if(any(colSums(m) == 1)){
      if(verbose){
        message(sprintf("Batch(es) %s only includes samples from one group of %s:",
                        paste0(which(colSums(m) == 1), collapse=','),
                        var))
        print(table(curr_batch_pid[,get(var)], curr_batch_pid[,batch]))
      }
      redo = T
    }
  }
  if(redo){
    if(verbose) message("Trying again to find more balanced batches...\n")
    curr_batch_pid = NULL
  }else{
    if(verbose) message("Success!")
  }
  return(list(batch_assignments = curr_batch_pid,
              success = !redo))
}


# make batches of roughly the same numbers of samples
# prioritize site and randgroupcode for randomization; check that age and sex are balanced by chance
make_batches_not_strict = function(nplates, curr_batch_pid, b, balance_vars){
  
  # while samples do not fit...
    # order samples by N and then group
    # for each row of curr_batch_pid
      # see if it fits in ANY batch
        # if yes, place and keep going 
        # if not, add plate; break
  
  curr_batch_pid[,group := paste0(codedsiteid, randomgroupcode)]
  curr_batch_pid = curr_batch_pid[order(N, group, decreasing = T)] 
  
  overflow = T
  while(overflow){
    
    # assign group, checking total size 
    batch_sizes = rep(0, nplates) 
    names(batch_sizes) = 1:nplates
    batch_iter = 1
    curr_batch_pid[,batch := NA_integer_]
    
    for (i in 1:nrow(curr_batch_pid)){
      
      # see if it fits in any batch
      if(!any(batch_sizes + curr_batch_pid[i, N] <= max_n_per_batch)){
        overflow = T
        if(verbose){
          message(sprintf("With this randomization, '%s' samples don't fit into %s batches. Trying with %s.",b, nplates,nplates+1))
        }
        nplates = nplates + 1
        overflow = T
        break
      }
      overflow = F
      
      # it fits in at least one batch. iterate over batches to put into first batch that it fits 
      which_fits = unname(which(batch_sizes + curr_batch_pid[i, N] <= max_n_per_batch))
      
      if(any(which_fits >= batch_iter)){
        # get first one greater than or equal to batch_iter
        which_fits = which_fits[which_fits >= batch_iter]
        which_fits = min(which_fits)
      }else{ 
        # if that doesn't work, start again at 1 
        which_fits[which_fits < batch_iter]
        which_fits = min(which_fits)
      }
      # add to batch
      curr_batch_pid[i, batch := which_fits] # assign batch
      batch_sizes[[which_fits]] = batch_sizes[[which_fits]] + curr_batch_pid[i, N] # increase size
      # increment batch
      batch_iter = batch_iter + 1
      if(batch_iter == nplates+1){
        batch_iter = 1
      }
    }
  }
  # at this point, all batches should be assigned
  batches = curr_batch_pid
  
  # check batch balance 
  check_batches = check_batch_balance(batches, balance_vars)
  if(!check_batches$success){
    # something, probably sex or age, was unbalanced
    # iteratively make batches until they are balanced
    warning("I didn't expect you to get here...")
    stop(sprintf("It looks like you want small batch sizes (N = %s). Please try re-running the script with the --strict-size flag for a batching method more suitable for small batch sizes.",
                 max_n_per_batch))
    batches = make_random_batches_not_strict(curr_batch_pid, b, max_n_per_batch, balance_vars)
  }
  
  return(batches)
}


id_optimal_batch_sizes = function(curr_batch_pid, max_n_per_batch){
  optimal_n_batches = ceiling(sum(curr_batch_pid[,N]) / max_n_per_batch)
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
  return(optimal_batch_sizes)
}


id_feasible_batch_sizes = function(curr_batch_pid, b, max_n_per_batch){
  
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
make_batches_strict = function(curr_batch_pid, b, max_n_per_batch, 
                               balance_vars = c('codedsiteid','randomgroupcode','sex_psca'),
                               verbose=T){
  
  feasible_batches = F
  
  n_samples_per_pid = curr_batch_pid[,N]
  names(n_samples_per_pid) = curr_batch_pid[,pid]
  n_samples_per_pid = as.list(n_samples_per_pid)
  
  # first try to find batches with ideal batch size 
  batch_sizes = id_optimal_batch_sizes(curr_batch_pid, max_n_per_batch)
  
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
      n_samples_per_pid_reordered = n_samples_per_pid_reordered[-1] 
      # put this in the first batch that it fits
      batch_to_place = as.numeric(names(remaining_room_per_batch)[remaining_room_per_batch>=set_to_place][1])
      if(is.na(batch_to_place)){
        break
      }
      new_sets[[batch_to_place]] = c(new_sets[[batch_to_place]], set_to_place)
      remaining_room_per_batch[[batch_to_place]] = remaining_room_per_batch[[batch_to_place]] - unlist(set_to_place)
    }
    # this worked? 
    if(all(remaining_room_per_batch==0)){
      if(verbose){
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
      check_batches = check_batch_balance(curr_batch_pid, balance_vars)
      balanced_batches = check_batches$success
      batch_assignments = check_batches$batch_assignments
      outer_loop = outer_loop + 1
      inner_loop = 1
    }else{
      inner_loop = inner_loop + 1
      if(inner_loop>1e6 & !feasible_batches){
        new_batch_sizes = id_feasible_batch_sizes(curr_batch_pid, b, max_n_per_batch)
        if(verbose){
          message(sprintf("With %s total samples and maximum %s samples per batch, the ideal number of samples per batch are as follows:\n%s\nHowever, after %s iterations, no combination of samples was found to fit these batch sizes. Trying again with the following batch sizes:\n%s",
                          sum(curr_batch_pid[,N]),
                          max_n_per_batch,
                          paste0(batch_sizes, collapse=', '),
                          inner_loop,
                          paste0(batch_sizes, collapse=', ')))
        }
        batch_sizes = new_batch_sizes
        feasible_batches = T
      }
      if(outer_loop > 200){
        if(verbose){
          writeLines(paste0(batch_sizes, collapse=', '))
        }
        stop(sprintf("With %s total samples, maximum %s samples per batch, and target batch sizes printed above, well-balanced batches were not found in %s candidate batches. You are currently requiring all batches to have samples from more than one group for all of the following variables:\n    %s\nTry removing the least important variable from this list using the --vars-to-balance flag and rerun the script.",
                     sum(curr_batch_pid[,N]),
                     max_n_per_batch,
                     outer_loop-1,
                     paste0(balance_vars, collapse=', ')))
      }
    }
  }
  return(batch_assignments)
}


make_random_batches_not_strict = function(curr_batch_pid, b, max_n_per_batch, 
                                          balance_vars = c('codedsiteid','randomgroupcode','sex_psca'),
                                          verbose=T){
  
  n_samples_per_pid = curr_batch_pid[,N]
  names(n_samples_per_pid) = curr_batch_pid[,pid]
  n_samples_per_pid = as.list(n_samples_per_pid)
  n_batches = ceiling(sum(curr_batch_pid[,N]) / max_n_per_batch)
  
  outer_loop = 1
  inner_loop = 1
  balanced_batches = F
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
      n_samples_per_pid_reordered = n_samples_per_pid_reordered[-1] 
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
      batch_iter = batch_iter + 1
      if(batch_iter > n_batches){
        batch_iter = 1
      }
    }
    # this worked? 
    if(length(n_samples_per_pid_reordered)==0){
      if(verbose){
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
      check_batches = check_batch_balance(curr_batch_pid, balance_vars)
      balanced_batches = check_batches$success
      batch_assignments = check_batches$batch_assignments
      outer_loop = outer_loop + 1
      inner_loop = 1
    }else{
      inner_loop = inner_loop + 1
      if(inner_loop>1e6){
        if(verbose){
          message("Can't find a combination of samples that fits in this many batches. Increasing the number of batches by 1.")
        }
        n_batches = n_batches + 1
      }
      if(outer_loop > 200){
        warning("I didn't expect you to get here either...")
        stop(sprintf("With %s total samples and up to %s samples per batch, well-balanced batches were not found in %s candidate batches. You are currently requiring all batches to have samples from more than one group for all of the following variables:\n    %s\nTry removing the least important variable from this list using the --vars-to-balance flag and rerun the script.",
                     sum(curr_batch_pid[,N]),
                     max_n_per_batch,
                     outer_loop-1,
                     paste0(balance_vars, collapse=', ')))
      }
    }
  }
  return(batch_assignments)
}


# CHECK FORMATS #########################################################################################################

# check formats
# if(length(shipments) == 0){
#   stop("Required argument '--shipment-manifest-excel' is empty. Please provide at least one path.")
# }
# if(length(apis) == 0){
#   stop("Required argument '--api-metadata-csv' is empty. Please provide at least one path.")
# }
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

# subset to existing samples
all_meta = all_meta[!(is.na(viallabel) | is.na(box))]

# SETUP #########################################################################################################

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
message(sprintf("These columns have 0 variance in the merged metadata: %s\n",
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
  message("'assay' is not in the column names of the merged API and shipment metadata. Batching will assume that all samples from a given tissue are for a single assay.\n")
  all_meta[,batching_group := sampletypecode]
}else{
  all_meta[,batching_group := paste0(assay, '_', sampletypecode)]
}

# BATCHING #########################################################################################################

for (b in unique(all_meta[,batching_group])){

  cat("_____________________________________________________________________________________________\n")
  cat(sprintf("--- BATCHING '%s' SAMPLES ---\n\n", b))
  
  curr_batch = unique(all_meta[batching_group == b])
  curr_batch_pid = unique(curr_batch[,.(codedsiteid, pid, randomgroupcode, sex_psca, calculatedage, older_than_40)])
  
  # how many samples per person?
  curr_batch_n = data.table(table(curr_batch[,pid]))
  colnames(curr_batch_n) = c('pid','N')

  curr_batch_pid[,pid := as.character(pid)]
  curr_batch_n[,pid := as.character(pid)]
  curr_batch_pid = merge(curr_batch_pid, curr_batch_n, by='pid')
  
  if(strict_size){
    batches = make_batches_strict(curr_batch_pid, b, max_n_per_batch, 
                                  balance_vars = balance_vars,
                                  verbose = verbose)
  }else{
    nplates = ceiling(nrow(curr_batch)/max_n_per_batch) 
    batches = make_batches_not_strict(nplates, curr_batch_pid, b, balance_vars) ## THIS IS GETTING STUCK
  }
  
  cat('Sample totals:\n')
  print(batches[,list(N_samples = sum(N)), by=batch])
  cat('\nN subj. per sex by batch:\n')
  print(table(batches[,batch], batches[,sex_psca]))
  cat('\nN subj. per site by batch:\n')
  print(table(batches[,batch], batches[,codedsiteid]))
  cat('\nN subj. per intervention by batch:\n')
  print(table(batches[,batch], batches[,randomgroupcode]))
  cat('\nN subj. per age group by batch (>40 yrs):\n')
  print(table(batches[,batch], batches[,older_than_40]))
  
  cat('\n\n')
  
  # write two versions to file 
  # batch characteristics 
  write.table(batches, file=sprintf("%s/precovid_%s-samples_UNBLINDED-batch-characteristics.csv", outdir, gsub(" ","-",b)), sep=',', col.names=T, row.names=F, quote=F)
  
  # current and new positions 
  curr_batch[,pid := as.character(pid)]
  all_info = merge(curr_batch, batches[,.(pid, batch)], by='pid')
  assert(nrow(all_info) == nrow(curr_batch))
  positions = all_info[,.(viallabel, barcode, sampletypecode, box, position, batch)]
  colnames(positions) = c('viallabel','barcode','sampletypecode','shipping_box','shipping_position','new_batch')
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

  write.table(positions, file=sprintf("%s/precovid_%s-samples_BLINDED-batch-assignments.csv", outdir, gsub(" ","-",b)), sep=',', col.names=T, row.names=F, quote=F)
}

warnings()
