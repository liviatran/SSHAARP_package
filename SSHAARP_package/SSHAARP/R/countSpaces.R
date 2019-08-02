#'Counts whitespaces in a string
#'
#'Provided string is split, and white spaces are counted in between split elements.
#'@param x A string.
#'
#'@note This function was obtained from Reddit and written by Josh Bredeweg, user jbraids1421.
#'
#'@return A numeric vector that identifies the number of spaces between each set of non-space characters in the input string.
#'
#'@export
#'
#'@examples countSpaces("The quick  brown   fox jumped    over  the lazy   dog.")
#
countSpaces <- function(x){
  counter <- 0
  coll <- numeric()
  vec <- strsplit(x," ")[[1]]
  for(i in 1:length(vec)){
    if (vec[i]==""){
      counter <- counter+1
    }
    else{
      if (counter!=0) coll <- c(coll,counter)
      counter <- 1
    }
  }
  coll
}
