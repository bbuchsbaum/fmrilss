#! /usr/bin/env Rscript



library(optparse)
library(fmrireg)

suppressPackageStartupMessages(library("optparse"))

option_list <- list(make_option(c("-b", "--bids_path"), type="character", help="the path to the BIDS folder containing the fMRI scans (default = '.')", default="."),
                    make_option(c("--bids_session"), type="character", help="the name of the bids session (optional)", default=""),
                    make_option(c("-d", "--deriv_folder"), type="character", help="the name of the subfolder containing the pre-processed scans (default = derivatives/fmriprep)", default="derivatives/fmriprep"),
                    make_option(c("--outdir"), type="character", help="the name of the output directory (default = betas)", default="betas"),
                    make_option(c("-t", "--task"), type="character", help="the name of the BIDS task"),
                    make_option(c("-l", "--duration"), type="numeric", default=1, help="the duration of each stimulus event (constant for all events; note: overrides 'duration' column in design file. default=1)"),  
                    make_option(c("-o", "--out"), type="character", default="betas", help="name of output file stem (default: betas)"),
                    make_option(c("-c", "--concatenate"), type="logical", default=FALSE, help="concatenate all beta estimates into one output file (default is FALSE)"),
                    make_option(c("-f", "--confounds"), type="character", help="a file containing a list of confound variables to be extracted from confounds.tsv files and added as nuisance regressors"),
                    make_option(c("-v", "--percvar"), type="numeric", help="percentage of confound variance to retain (default is 95)", default=95),
                    make_option(c("-e", "--est"), type="character", default="lss", help="method for beta estimation ('lss' or 'lsa'; default is 'lss'"),
                    #make_option(c("--censor"), type="character", help="list of censor files, one per run each on a separate line"),
                    make_option(c("--milliseconds"), type="logical",  action="store_true", default=FALSE, help="onsets are in milliseconds"),
                    make_option(c("--onsets"), type="character", help="name of onset column (default is 'onsets'", default="onsets"),
                    #make_option(c("--fixed_onsets"), type="character", help="name of onset column of fixed event set"),
                    #make_option(c("--fixed_factor"), type="character", help="name of column of fixed variable"),
                    #make_option(c("-r", "--runs"), type="character", help="optional name of file containing runs to subselect from design file"),
                    #make_option(c("--runcolumn"), type="character", help="name of column in design file indicating the run/block number", default="runs"),
                    make_option(c("-p", "--polort"), type="numeric", default=3, help="number of polynomial regressors -- passed to as 'polort' 3dDeconvolve (default = 3)"),
                    make_option(c("-m", "--mask"), type="character", help="name of binary image mask"),
                    make_option(c("-s", "--subid"), type="character", help="the subject id"),
                    make_option(c("-x", "--space"), type="character", help="the name of the pre-processed bold space", default="MNI152NLin2009cAsym_preproc"))
                   
oparser <- OptionParser(usage = "estimate_betas2.R [options]", option_list=option_list)
opt <- parse_args(oparser, positional_arguments=TRUE)
args <- opt$options



## load in confounds variable list
cvars <- if (!is.null(args$confounds)) {
             ret <- scan(args$confounds, "")
             message("confound variables are: ", paste(ret, collapse=" "))
             ret
         }

if (args$bids_path == ".") {
    args$bids_path <- getwd()
}

## get bids-fmriprep information

session <- if (args$bids_session != "" && !is.null(args$bids_session)) {
   args$bids_session
} else {
    ""
}
               
print(session)
bids <- bidser::bids_project(args$bids_path, fmriprep=TRUE)
#bids <- fmrireg:::bids_source(args$bids_path, args$deriv_folder, id=args$subid, bold_space=args$space, task=args$task, session=session, confound_vars=cvars)

                              
if (!is.null(cvars)) {
    ## TODO deal with sessions
    ## read in counds and retain 'perc_var' as PCs
    cdat <- bidser:::read_confounds(bids, perc_var=args$percvar, cvars=cvars, subid=args$subid, task=args$task, nest=TRUE)
    #nuisanceCovars <- split(cdat$data, cdat$data$run)
    nuisanceCovars <- lapply(cdat$data, as.matrix)
    nuisanceCovars <- lapply(nuisanceCovars, function(x) {
        keep <- !apply(x, 2, function(z) all(is.na(z)))
        x[,keep, drop=FALSE]
        })
} else {
  nuisanceCovars = NULL
}

## construct vector of scan filenames
scans <- paste0(args$bids_path, "/", bidser:::preproc_scans(bids, subid=args$subid, task=args$task))
#scans <- paste0(bids$func_path, "/", bids$preproc_scans)
fexists <- sapply(scans, file.exists)

if (!all(fexists)) {
  stop(paste("could not find input file(s):", paste(scans[!fexists], collapse=" ")))
}

cdir <- getwd()
outpath <- file.path(dirname(dirname(scans[1])), args$outdir)

if (!file.exists(outpath)) {
    system(paste("mkdir ", outpath))
}

message("setting current working directory to:", outpath)
setwd(outpath)

## name of onset column in events files
onsetColumn <- args$onsets

if (args$milliseconds) {
    message("onsets are in milliseconds")
    #design[[onsetColumn]] <- design[[onsetColumn]]/1000
}

#censorFiles <- if (!is.null(args$censor)) {
#    cfiles <- scan(args$censor,"")
#    if (length(cfiles) != length(runLevs)) {
#        stop("must have 1 censor file for each scanning run")
#    }
   
#    res <- lapply(cfiles, function(fname) read.table(fname, header=TRUE)[,1])
#    cnames <- paste0("censorVec_", runLevs, ".1D")
#    
#    for (i in 1:length(runLevs)) {
#        censorVec <- res[[i]]
#        write(censorVec,cnames[i], ncolumns=1)
#    }
#
#    cnames
#}

censorFiles <- NULL

#if (!is.null(design$duration)) {
#  message("found duration column in design file: will use for event duration vector")
#  duration <- design$duration
#} else {
#  duration <- rep(as.numeric(args$duration), nrow(design))
#  message("using constant event duration of: ", duration[1])
#}

duration <- as.numeric(args$duration)
message("using constant event duration of: ", duration)

#print(args)
#print(bids$preproc_scans)


## TODO generate reasonable output names for new fmriprep ....
outfiles <- basename(gsub("desc-preproc", "desc-betas", scans))
#outfiles <- past0("sub-", args$subid, "_task-", args$task, "_run-", 


message("outfiles: ", paste(outfiles, collapse= " "))

if (!is.null(args$mask)) {
    if (is.null(session) || session == "") {
        args$mask <- file.path(args$bids_path, args$deriv_folder, paste0("sub-", args$subid), "func", paste0("sub-", args$subid, "_", args$mask))
    } else {
        m <- file.path(args$bids_path, args$deriv_folder, paste0("sub-", args$subid), paste0("ses-", session), "func", paste0("sub-", args$subid, "_ses-", session, "_", args$mask))

        if (!file.exists(m)) {
            m <- file.path(args$bids_path, args$deriv_folder, paste0("sub-", args$subid), paste0("ses-", session), "func", paste0("sub-", args$subid, "_", args$mask))
        }

        args$mask <- m
        
    }
    message("mask is: ", args$mask)
}

genStimFileArgs <- function(fnames) {
  
  paste(lapply(1:length(fnames), function(k) {
    paste("-stim_base", (k+1), "-stim_file", (k+1), fnames[k], "-stim_label", (k+1), paste0("nuisance_", k))
  }), collapse = " ")
}

#genFixedFileArgs <- function(fnames) {
#}
  

message("starting estimation")


#evfiles <- bidser:::event_files(bids, subid=args$subid, task=args$task)
evfiles <- bidser:::search_files(bids, subid=args$subid, task=args$task, kind="events", full_path=TRUE)

for (run in seq_along(scans)) {
    message("run number: ", run)
    message("scan: ", scans[run])
    
  nuisnames <- if (!is.null(nuisanceCovars)) {
    nuis <- nuisanceCovars[[run]]
    sapply(1:ncol(nuis), function(j) {
      oname <- paste0("nuisance_", args$task, "run_", run, "_reg#", j, ".1D")
      write(nuis[,j], file=oname, ncolumns=1)
      oname
    })
  } else {
    NULL
  }

  nstim <- if (!is.null(nuisanceCovars)) {
    length(nuisnames) + 1
  } else {
    1
  }

  design <- read.table(evfiles[run], header=TRUE, na.strings="n/a")

  if (is.null(design[[onsetColumn]])) {
    stop(paste("design file must have column named:", onsetColumn))
  }

  onsets <- design[[onsetColumn]]

  if (args$milliseconds) {
      onsets <- onsets/1000
  }

  rundurs <- rep(duration, nrow(design))
  onsetMat <- cbind(onsets, rundurs)
  onsetStr <- sapply(1:nrow(onsetMat), function(i) paste0(onsetMat[i,1], ":", onsetMat[i,2]))
  

  onsname <- paste0("onsets_", args$task, "run", run, ".1D")

  message(onsetStr)
  write(onsetStr, file=onsname, ncolumns=length(onsets))

  dmat <- paste0("R.", args$task, "xmat.", run, ".1D")

  if (args$est == "lsa") {

    cmd1 <- paste("3dDeconvolve -input",
                scans[run],
                if (!is.null(args$mask)) paste("-mask", args$mask, collapse=" ") else "-automask",
                "-polort", args$polort, "-x1D_uncensored", dmat,              
                "-num_stimts", nstim,
                "-stim_times_IM 1", onsname, paste0("dmBLOCK"),
                "-stim_label 1 event",
                "-nobout",
                "-bucket", outfiles[run],
                  "-nofullf_atall",
                if (!is.null(censorFiles)) paste("-censor", censorFiles[run]) else "",
                if (!is.null(nuisnames)) genStimFileArgs(nuisnames) else "")
      print(paste("RUNNING:", cmd1))
      
      system(cmd1)
      cmd2 <- paste("3dTcat -prefix tmp.nii.gz", outfiles[run])
      system(cmd2)
      system("rm tmp.nii.gz")
      
  } else {

      cmd1 <- paste("3dDeconvolve -input",
                scans[run],
                "-polort", args$polort, "-x1D_uncensored", dmat,              
                "-x1D_stop -num_stimts", nstim,
                "-stim_times_IM 1", onsname, paste0("dmBLOCK"),
                "-stim_label 1 event",
                if (!is.null(censorFiles)) paste("-censor", censorFiles[run]) else "",
                if (!is.null(nuisnames)) genStimFileArgs(nuisnames) else "")
                
   

      cmd2 <- paste("3dLSS -verb", "-prefix", outfiles[run], "-input", scans[run],
                    if (!is.null(args$mask)) paste("-mask", args$mask, collapse=" ") else "-automask",
                    "-matrix", dmat)
      print(paste("RUNNING:", cmd1))
      system(cmd1)
      print(paste("RUNNING", cmd2))
      system(cmd2)
  }
}




if (args$concatenate) {
    conc_out <- paste0(args$out, "_all.nii.gz")
    cmd <- paste("3dTcat -prefix", conc_out, paste(outfiles, collapse=" "))
    system(cmd)
}



#if (!is.null(rootdir)) {
#  system(paste("rm -R", rootdir))
#}
