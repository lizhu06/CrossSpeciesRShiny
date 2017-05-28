upperlowerMatch <- function(upperset,lowerset){
  #upper/lower case matching
  match_results <- match(toupper(lowerset),upperset)
  lowerset_match <- lowerset[which(!is.na(match_results))]
  upperset_match <- upperset[match_results[!is.na(match_results)]]
  
  upperset_unmatch <- setdiff(upperset,upperset_match)
  lowerset_unmatch <- setdiff(lowerset,lowerset_match)
  
  out <- list(upperset_match,lowerset_match,
              upperset_unmatch,lowerset_unmatch)
}  