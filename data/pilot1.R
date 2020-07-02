data.path=file.path("data", "raw", "pilot1")

data.files=list.files(path=data.path, pattern=".*_pilot1_.*\\.csv", full.names = T )

pilot1 <- do.call(rbind, 
                  lapply(data.files, function(fname){
                    read.table(fname, sep=",", header=T, comment.char = "#", stringsAsFactors = F);
                  }))

# fix subj 11 who used left/right for exactly 1 response (??)
#pilot1 %>% filter(stimulus==0, subj==11) %>% do(broom::tidy(table(.$response)))
pilot1 <- within(pilot1, {
  response[stimulus==0 & response=="left"]="lalt";
  response[stimulus==0 & response=="right"]="rctrl";
})
