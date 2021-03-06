# 2020-ftrsgt-paper

This repository hosts code associated with the paper

The interplay between cognitive control, behavioral variability and mind wandering: Insights from a HD-tDCS study https://psyarxiv.com/d9ngb/


## Format

The data is prepared using the [ProjectTemplate](http://projecttemplate.net/). In this format, all raw data is stored in `data` and automatically loaded into R-datasets when running

~~~{R}
library(ProjectTemplate)
load.project()
~~~

in an open R-session. It is possible that you will have to run `migrate.project()` the first time you load the datasets like so:

~~~{R}
library(ProjectTemplate)
migrate.project()
load.project()
~~~

The available datasets from the three studies are `pilot1`, `pilot2` and `stimulation`:

~~~{R}
> str(pilot1)
'data.frame':	99488 obs. of  6 variables:
 $ subj     : int  1 1 1 1 1 1 1 1 1 1 ...
 $ sessionix: int  1 1 1 1 1 1 1 1 1 1 ...
 $ ISI      : num  1.25 1.25 1.25 1.25 1.25 1.25 1.25 1.25 1.25 1.25 ...
 $ time     : num  0.00156 0.37991 0.92388 1.25036 2.50033 ...
 $ stimulus : int  1 0 0 1 1 0 1 0 1 0 ...
 $ response : chr  "" "lalt" "rctrl" "" ...
> str(pilot2)
'data.frame':	67374 obs. of  5 variables:
 $ subj    : int  1 1 1 1 1 1 1 1 1 1 ...
 $ trial   : int  0 1 1 2 2 3 3 4 4 5 ...
 $ time    : num  0.0013 0.7506 0.9522 1.5005 1.5922 ...
 $ stimulus: chr  "stimulus" "stimulus" "tap" "stimulus" ...
 $ response: chr  "" "" "rctrl" "" ...
> str(stimulation)
'data.frame':	191266 obs. of  5 variables:
 $ subj    : int  1 1 1 1 1 1 1 1 1 1 ...
 $ trial   : int  0 0 1 1 2 2 3 3 4 4 ...
 $ time    : num  0.00014 0.54114 0.75017 1.03724 1.50024 ...
 $ stimulus: chr  "stimulus" "tap" "stimulus" "tap" ...
 $ response: chr  "" "rctrl" "" "lctrl" ...
 ~~~

 A pre-processed dataset can be created by running

 ~~~{R}
get.nback(d = stimulation)
 ~~~