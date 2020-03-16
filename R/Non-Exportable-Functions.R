## ---- Non-Exportable-Functions.R ---- ##
#                                        #
#     .get.exepath                       #
#     .test.exefile                      #
#     msa                                #
#                                        #
## ------------------------------------ ##

## --------------------------------------------------------------- ##
#                      .get.exepath(prg)                            #
## --------------------------------------------------------------- ##
#' Find Full Paths to Executables
#' @description Finds the path to an executable.
#' @param prg name of the executable.
#' @return Returns the absolute path.
#' @keywords internal.

.get.exepath <- function(prg) {

  paths <- list(

    pymol = list(

      Linux = c("/usr/bin/pymol",
                "/usr/local/bin/pymol"),

      Darwin = c("/Applications/MacPyMOL.app/Contents/MacOS/MacPyMOL",
                 "/Applications/MacPyMOLX11Hybrid.app/Contents/MacOS/MacPyMOL",
                 "/usr/bin/pymol",
                 "/usr/local/bin/pymol"),

      Windows = c("C:/python27/PyMOL/pymol.exe",
                  "C:/Program Files/PyMOL/PyMOL/PymolWin.exe",
                  "C:/Program Files/PyMOL/PymolWin.exe"),

      ver = "-cq"
    ),


    muscle = list(

      Linux = c("/usr/bin/muscle",
                "/usr/local/bin/muscle"),

      Darwin = c("/usr/bin/muscle",
                 "/usr/local/bin/muscle"),

      Windows = c("C:/Program Files/muscle.exe",
                  "C:/Program Files/muscle3.8.31_i86win32.exe",
                  "C:/Program Files/muscle/muscle.exe",
                  "C:/Program Files/Muscle/muscle.exe",
                  "C:/Program Files/seaview/muscle.exe",
                  "C:/Program Files/seaview4/muscle.exe"),

      ver = "-version"
    ),

    foldx = list(

      Linux = c("/usr/bin/foldx",
                "/usr/local/bin/foldx"),

      Darwin = c("/usr/bin/foldx",
                 "/usr/local/bin/foldx"),

      Windows = c("C:/Program Files/foldx.exe",
                  "C:/Program Files/foldx/foldx.exe",
                  "C:/Program Files/FoldX/foldx.exe"),

      ver = "--version"
    ),

    dssp = list(

      Linux = c("/usr/bin/dssp",
                "/usr/local/bin/dssp"),

      Darwin = c("/usr/bin/dssp",
                 "/usr/local/bin/dssp",
                 "/usr/bin/mkdssp",
                 "/usr/local/bin/mkdssp",
                 "/anaconda3/bin/dssp",
                 "/anaconda3/bin/mkdssp"),

      Windows = c("C:/Program Files/dssp.exe",
                  "C:/Program Files/dssp-2.0.4-win32.exe",
                  "C:/Program Files/dssp/dssp.exe",
                  "C:/Program Files/Dssp/dssp.exe"),

      ver = "--version"
    )

  )

  ## user provided full path
  if(file.exists(prg) & !dir.exists(prg)) {
    return(prg)
  }

  ## try to automatically determine path
  exefile <- Sys.which(prg)

  if(nchar(exefile) == 0) {

    if(prg %in% c("pymol", "muscle", "clustalo", "dssp")) {
      ## determine os
      os1 <- Sys.info()["sysname"]

      ## use guess-paths defined above
      exefiles <- paths[[prg]][[os1]]
      fe <- file.exists(exefiles)

      if(any(fe)) {
        exefile <- exefiles[which(fe)[1]]
      }
      else {
        exefile <- NULL
      }
    }
    else {
      exefile <- NULL
    }
  }

  if(is.null(exefile)) {
    stop(paste0("could not determine path to '", prg, "'"))
  }
  return(exefile)
}

## --------------------------------------------------------------- ##
#                      .test.exefile(exefile)                       #
## --------------------------------------------------------------- ##
.test.exefile <- function(exefile) {
  prg <- tolower(basename(exefile))

  if(grepl("muscle", prg)) {
    ver <- "-version"
  }
  if(grepl("pymol", prg)) {
    ver <- "-cq"
  }
  if(grepl("foldx", prg)) {
    ver <- "--version"
  }
  if(grepl("dssp", prg)) {
    ver <- "--version"
  }


  os1 <- Sys.info()["sysname"]
  if (os1 == "Windows") {
    success <- shell(paste(shQuote(exefile), ver))
  }
  else {
    success <- system(paste(exefile, ver),
                      ignore.stderr = TRUE, ignore.stdout = TRUE)
  }

  if(!(success %in% c(0,1))) {
    return(FALSE)
  }
  else {
    return(TRUE)
  }
}

## ----------------------------------------------------------------- ##
#   msa <- function(sequences, ids = names(sequences), sfile = FALSE) #                #
## ----------------------------------------------------------------- ##
#' Multiple Sequence Alignment
#' @description Aligns multiple protein sequences
#' @param sequences vector containing the sequences.
#' @param ids vector containing the sequences' ids.
#' @param sfile path to the file where the fasta alignmet should be saved, if any.
#' @return Returns a list of four elements. The first one ($seq) provides the sequences analyzed, the second element ($ids) retuns the identifiers, ther third element ($aln) privides the alignment in fasta format and the fourth element ($ali) gives the alignment in matricial format.
#' @keywords internal.

msa <- function(sequences, ids = names(sequences), sfile = FALSE){

  if (requireNamespace('Biostrings', quietly = TRUE)){
    seq <- Biostrings::AAStringSet(sequences)
  } else {
    stop("You must install the package Biostrings in order to use this function")
  }

  if (requireNamespace('muscle', quietly = TRUE)){
    aln1 <- muscle::muscle(seq)
  } else {
    stop("You must install the package muscle in order to use this function")
  }

  aln <- list()
  aln$seq <- sequences
  aln$ids <- ids
  aln$aln <- as.character(aln1)
  l <- sapply(aln$aln, function (x) strsplit(x, split = ""))
  aln$ali <- matrix(unlist(l), nrow = length(sequences), byrow = TRUE)

  if (sfile != FALSE){
    for (i in 1:length(aln$aln)){
      t <- paste(">", aln$ids[i], sep = "")
      cat(t, file = sfile, append = TRUE)
      tt <- paste("\n", aln$aln[i], "\n", sep = "" )
      cat(tt, file = sfile, append = TRUE)
    }
  }
  return(aln)
}
