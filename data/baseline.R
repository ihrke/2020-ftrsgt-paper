data.path=file.path("data", "raw", "hdtdcs")

data.files=list.files(path=data.path, pattern=".*_baseline_.*\\.csv", full.names = T )

baseline <- do.call(rbind,
                     lapply(data.files, function(fname){
                       read.table(fname, sep=",", header=T, comment.char = "#", stringsAsFactors = F);
                     }))



