data.path=file.path("data", "raw", "pilot1")

data.files=list.files(path=data.path, pattern=".*_RNGT_.*\\.csv", full.names = T )

pilot1.rngt <- do.call(rbind, 
                 lapply(data.files, function(fname){
                   read.table(fname, sep=",", header=T, comment.char = "#", stringsAsFactors = F);
                 }))

