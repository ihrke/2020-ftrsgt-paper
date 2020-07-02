data.path=file.path("data", "raw", "pilot2", "behavior")

data.files=list.files(path=data.path, pattern=".*_pilot4_.*\\.csv", full.names = T )

pilot2 <- do.call(rbind, 
                  lapply(data.files, function(fname){
                    read.table(fname, sep=",", header=T, comment.char = "#", stringsAsFactors = F);
                  }))

# remove two tap-responses that where "return" and "right" for some reason (instead of lalt and rctrl)
pilot2=within(pilot2, {
  response[response=="return"]="lalt"
  response[response=="right"]="rctrl"
})



d=get.nback(pilot2, nback=20, which.apen.m = 1, on.task.crit = 2)
d %<>% select(-part,-condition)
