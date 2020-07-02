data.path=file.path("data", "raw", "hdtdcs")

data.files=list.files(path=data.path, pattern=".*_stimulation_.*\\.csv", full.names = T )

stimulation <- do.call(rbind,
                     lapply(data.files, function(fname){
                       read.table(fname, sep=",", header=T, comment.char = "#", stringsAsFactors = F);
                     }))

