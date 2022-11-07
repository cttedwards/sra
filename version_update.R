versionUpdate <- function(version, update = "0.0.0") {
  
  version <- as.integer(unlist(strsplit(version, split = ".", fixed = TRUE)))
  update  <- as.integer(unlist(strsplit(update,  split = ".", fixed = TRUE)))
  
  return(paste(version + update, collapse = "."))
}

# Get Version and date
VERSION <- scan('DESCRIPTION',what = character(),skip = 3,nlines = 1)[2]
VERSION <- versionUpdate(VERSION, "0.0.1")
DATE    <- Sys.Date()
TIME    <- Sys.time()

# update DESCRIPTION
DESCRIPTION    <- readLines('DESCRIPTION')
DESCRIPTION[4] <- paste('Version:', VERSION)
DESCRIPTION[5] <- paste('Date:', DATE)
writeLines(DESCRIPTION, 'DESCRIPTION')
rm(DESCRIPTION)

# Write .onAttach
filename <- "R/zzz.R"
cat(".onAttach <- function(libname, pkgname)\n", file = filename)
cat("{\n", file = filename, append = TRUE)
cat(paste("    packageStartupMessage(\"sra v.", VERSION, " (", TIME, ");\\nInputs compiled using: sraInputsBio v.", packageVersion("sraInputsBio"), "; sraInputs v.", packageVersion("sraInputs"), "\")\n", sep = ""), file = filename, append = TRUE)
cat("}\n", file = filename, append = TRUE)
rm(filename)
