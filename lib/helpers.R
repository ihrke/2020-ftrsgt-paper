if("MASS" %in% loadedNamespaces()){
  select <- dplyr::select
}


## format p-value
pformat <- function(val) { sub("^(-?)0.", "\\1.", sprintf("%.2f", val)) }


roll.apen <- function(x, nback=15, maxM=3){
  res = matrix(NA, nrow=length(x), ncol=maxM+1)
  for(i in nback:length(x)){
    v=x[(i-nback+1):i]
    if(any(is.na(v))){
      res[i,]=NA
    } else {
      res[i,] = apen_int(v, maxM)
    }
  }
  res
}
roll.sd <- function(x, nback=15){
  res = vector(length=length(x))*NA
  for(i in nback:length(x)){
    res[i] = sd(x[(i-nback):i], na.rm=T)
  }
  res
}
## modified from: http://stackoverflow.com/a/11813377/972791
save.append <- function(..., list = character(), file) {
  if(file.exists(file)){
    previous  <- load(file)
  } else {
    previous=c()
  }
  var.names <- c(list, as.character(substitute(list(...)))[-1L])
  for (var in var.names) assign(var, get(var, envir = parent.frame()))
  save(list = unique(c(previous, var.names)), file = file)
}

## pretty print with NAs
print.df.na <- function(df, na.label='.'){
  df <- data.frame(df)
  .tmp <- format(df)
  .tmp[is.na(df)] <- na.label
  return(.tmp)
}

install.common.packages <- function(){
  needed_pkgs=c('mvtnorm', 'ggplot2', 'dplyr', 'effects', 'pwr', 'arm', 'psych', 'car', "QuantPsyc",
                "rjags", "reshape2", "tidyr", 'mvtnorm', 'truncnorm', 'Rcpp', 'ProjectTemplate',
                'devtools', 'rstan')

  # Install CRAN packages (if not already installed)
  .inst <- needed_pkgs %in% installed.packages()
  if(length(needed_pkgs[!.inst]) > 0) install.packages(needed_pkgs[!.inst])

  # Load packages into session 
  #.tmp=lapply(needed_pkgs, require, character.only=TRUE)
}


printf <- function(s, ...){
  cat(sprintf(s, ...))
}

sem <- function(x){
  (sd(x)/sqrt(length(x)))
}

softmax <- function(x){
  if(is.null(dim(x))){
    dim(x) <- c(1,length(x))
  }
  r<-exp(x)/rowSums(exp(x))
  if(dim(x)[1]==1){
    r<-as.vector(r)
  }
  r
}

## find index in x that is closest to val
closest <- function(x, val){
  which.min(abs(x-val))
}
