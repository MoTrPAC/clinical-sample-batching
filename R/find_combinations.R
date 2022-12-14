#' Combination sum
#' 
#' Return a list of all possible combinations of \code{numbers} that sum to \code{target}. 
#' Each value in \code{numbers} can be used unlimited times. 
#'
#' @param numbers vector of possible integers
#' @param target integer, target sum
#'
#' @return list of vectors
#' @export 
#'
#' @examples
#' combination_sum(c(2,3,5), 8)
#' combination_sum(c(1,2,3,5), 8)
#' 
#' @source <https://leetcode.com/problems/combination-sum/solutions/16506/8-line-python-solution-dynamic-programming-beats-86-77/>
combination_sum = function(numbers, target){
  numbers = numbers[order(numbers, decreasing = FALSE)]
  
  # make a list of lists of vectors 
  dp = list()
  for(i in 1:target){
    dp[[i]] = list()
  }
  
  for(i in 1:target){
    for(n in numbers){
      if(n > i){
        break
      }
      if(n == i){
        dp[[i]] = append(dp[[i]], list(n))
      }else{
        if(length(dp[[i-n]]) > 0){
          for(j in 1:length(dp[[i-n]])){
            currnums = dp[[i-n]][[j]]
            if(length(currnums) == 1){
              lastnum = currnums
            }else{
              lastnum = currnums[length(currnums)]
            }
            if(n >= lastnum){
              dp[[i]] = append(dp[[i]], list(c(dp[[i-n]][[j]], n)))
            }
          }
        }
      }
    }
  }
  return(dp[[target]])
}

target_combinations = function(possible_combinations, curr_batch_pid){
  counter = list()
}


