#!/usr/bin/env Rscript

## ========================================== ## 
## LOAD PACKAGES                             ##
## ========================================== ## 
library(optparse)
library(dplyr)


## ========================================== ## 
## PARSE COMMAND-LINE PARAMETERS              ##
## ========================================== ## 

option_list <- list(
  make_option(c("-i", "--sample_sheet"  ), type="character" , default=NULL   , metavar="path"   , help="Sample Sheet, colData, where samples are rows and columns are sample annotations."    ), 
  make_option(c("-b", "--iem_format"    ), type="logical"   , default=TRUE   , metavar="boolean", help="If Sample sheet is in illumina IEM format."                                            ),
  make_option(c("-a", "--require_adapter"), type="logical"  , default=TRUE  , metavar="boolean", help="If to check for Illumina adapter. Will only check if there is a character vector here or not"),
  make_option(c("-d", "--allow_dup_snames"    ), type="logical"   , default=FALSE  , metavar="boolean", help="If to check for duplicate sample names."                      ),
  make_option(c("-u", "--allow_num_prefix"    ), type="logical"   , default=FALSE  , metavar="boolean", help="If to allow sampple names to start with numeric character. This may complicate downstream R analysis that do not allow column names to start with numerics."                      ),
  make_option(c("-o", "--output_ssheet"    ), type="logical"   , default=TRUE  , metavar="boolean", help="If to output a sample sheet where forbidden characters are changed to underscore-"    ),
  make_option(c("-n", "--output_nfcore"    ), type="logical"   , default=TRUE  , metavar="boolean", help="If to output nfcore style sample-sheet-"                          ),
  make_option(c("-p", "--paired"        ), type="logical"   , default=TRUE   , metavar="boolean", help="If pariend end or not."                                            ),
  make_option(c("-s", "--strandness"    ), type="character"   , default="reverse"   , metavar="string", help="If pariend end or not."                                       ),
  make_option(c("-m", "--bam_suffix"    ), type="character"   , default="Aligned.sortedByCoord.out.bam"   , metavar="string", help="Suffix to use togehter with Sample_Names when generating matching BAM file names in sample sheet. The BAM suffix is currently defined in  STAR"),
  make_option(c("-c", "--count_mat_siuffix"    ), type="character"   , default="Aligned.sortedByCoord.out.bam"   , metavar="string", help="Suffix to use togehter with Sample_Names when generating matching Count matrix sample id. This suffix is currently defined in by STAR and featureCounts")
)


opt_parser <- OptionParser(option_list=option_list, add_help_option=FALSE)
opt        <- parse_args(opt_parser)

if (is.null(opt$sample_sheet)){
  print_help(opt_parser)
  stop("Please provide a sample_sheet file.", call.=FALSE)
}

# check strandness
strandness.types <- c("reverse","forward","unstranded")
if(!opt$strandness %in% strandness.types) stop("--strandness must be one of: ", paste(strandness.types, collapse=", "))

cat("\n ... CHECKING SAMPLE SHEET\n\n")


## ========================================== ## 
## PROCESS COUNTS FILE  - IEM STYLE
## ========================================== ## 
## comment table is read and hashes are commented, hashes used when featureCounts produce output

#ssheet <- read.delim()
# opt$sample_sheet <- "/Users/david/CTG_projects/bulkRNA/CTG_2020_177/samplesheet_2020_177.csv"
all.lines <- scan(file = opt$sample_sheet, what = "character", sep="\n", nmax = 250)

## check nlines of document
if(length(all.lines)>=250) stop("Number of rows in sample sheet exceeds what is read by the scan function ('nmax' set to 250)")

## check file delimiter, column separators
# u <- grep(opt$sep, all.lines[1])
# if(!length(u)) stop("Cant find separators in string, check if correct file delimiter format is used using --sep")

## check IEM sections
if(opt$iem_format){
  iem.headers <- list(Header="[[]Header[]]", Reads="[[]Reads[]]",Settings="[[]Settings[]]",Data="[[]Data[]]")
  iem.index <- sapply(iem.headers, function(x) grep(x, all.lines))
  if(!any(unlist(lapply(iem.index, length)==1))){
    stop("Sample Sheet Not IEM format - check Illumina IEM sections, must include - [Header], [Reads], [Settings], [Data] ...\n ... or set --iem_format to FALSE")
  }
  iem.index <- sapply(iem.headers, function(x) grep(x, all.lines))
  
  ## Check Adapter if adapter sequence slot is NA or not
  # Note, will not check if different Adapter reads for forward or reverse
  if(opt$require_adapter){
    settings.section <- read.delim(file=opt$sample_sheet, header = F, sep = ",", 
                                   skip = iem.index["Settings"], nrows = iem.index["Data"]-iem.index["Settings"]-1)
    u <- grep("Adapter", settings.section[,1])
    if(!length(u)) stop("Cant find 'Adapter' row in [Settings] section ")
    if(is.na(settings.section[u[1],1])) stop("No Adapter sequence detected. Define adapter or set --require_adapter FALSE")
 }
}else{
  iem.index = c(Data=0)
}
  
  # read data sheet as data frame (from below Data section header)
  data_table <- read.delim(file=opt$sample_sheet, header = TRUE, sep = ",", 
                             skip = iem.index["Data"])
  
  
  # Check Data section. 
  ## 1. Columns must include Sample_ID,Sample_Name,Sample_Project columns
  required.columns <- c("Sample_ID","Sample_Name","Sample_Project")
  if(!all (required.columns %in% colnames(data_table))){
    stop("not all required header names are present")
  }



# 2a. Check for non permitted special characters
## The Description column, if present, is treated a bit more gentle. Only type_a specials are disallowed here amd spaces, dashes, parenthese, dots are allowed
specials <- list(
    type_a=c("[*]","[/]","[+]","[']","[`]","[?]","[=]"),
    type_b=c("[:]","[-]","[(]","[)]","[.]"," ")
    )

# matrix(data_table)
special.flag <- F
data_table.out <- data_table

cat("\n\n ... ** CHECKING DATA COLUMNS FOR SPECIAL CHARACTERS ** ")
special.flag <- F # used to keep track if special characters are present in any column

for(i in 1:ncol(data_table.out)){
  # test_specials <- unlist(specials)
  # if(colnames(data_table.out)[i]=="Description") test_specials <- specials$type_a
  ## the required.columns will be tested for the most strict character list.
  if(colnames(data_table.out)[i] %in% required.columns){ 
    test_specials <-  unlist(specials)
  }else{
    test_specials <- specials$type_a
    }
  
  
  for(ii in 1:length(test_specials)){
    u<-grepl(test_specials[ii], data_table.out[,i]) 
    if(any(u)){ 
      special.flag <- T
      cat("\n ... ... Warning: column  '", colnames(data_table.out)[i], "'  contains: ", ifelse(test_specials[ii]==" ", "white space", test_specials[ii]))  
      ## replace these with underscore in out-file 
      data_table.out[,i] <- gsub(test_specials[ii], "_", data_table.out[,i])
      }
  }
}

## Check if Sampl_ID contain duplicates
u <- duplicated(data_table$Sample_ID)
if(any(u)) stop("Duplicated Sample_ID present: ", data_table$Sample_ID[u])
cat("\n\n ... OK, no duplicate Sample_IDs found.")



## Check if Sampl_Names contain duplicates
if(!opt$allow_dup_snames){
  u <- duplicated(data_table$Sample_Name)
  if(any(u)) stop("Duplicated Sample_Name present: ", data_table$Sample_Name[u])
  cat("\n\n ... OK, no duplicate Sample_Names found.")
}

# check if any Sample_Name starts with numeric character
if(!opt$allow_num_prefix){
  chr_vec <- sapply(data_table$Sample_Name, function(x){
    y <- substr(x, 1, 1)
    !is.na(suppressWarnings(as.numeric(y)))
  })
if(any(chr_vec)) stop("There are Sample_Name(s) starting with numeric characters: ", paste(my_samples[chr_vec], collapse = "; "))
}

# If to save updated sample sheet to use in next run/itteratve process
if(opt$output_ssheet & special.flag){
  out.file.name <- paste(opt$sample_sheet, "_iem_corrected.csv", sep="")
  if(grepl("[.]csv", opt$sample_sheet)) out.file.name <- gsub(pattern = "[.]csv", replacement = "_iem_corrected.csv", x = opt$sample_sheet)
  # if(file.exists(out.file.name)) stop("Out file already exists: ",out.file.name)
  cat("\n\n ... writing corrected file to output: ",out.file.name)
  
  cat(all.lines[1:iem.index["Data"]], file = out.file.name, sep = "\n", append = F)
  cat(colnames(data_table.out), file = out.file.name, sep = ",", append = T)
  cat("\n", file = out.file.name, sep = ",", append = T)
  write.table(x = data_table.out, file = out.file.name, sep = ",", append = T, quote = F, row.names = F, col.names = F, na="", )
}

## stop with error if non-allowed special characters  were found
if(special.flag) stop("Sample sheet contains special characters. Fix this !!")



cat("\n ... ### --------------------### " )
cat("\n ... ### Sample Sheet OK  " )
cat("\n ... ### ------------------- ### \n\n" )


if(opt$output_nfcore){
  cat("\n\n\n ... \n ...  ......\n" )
  cat("\n ... ### -------------------------------------"      )
  cat("\n ... ###   GENERATE NF-CORE SAMPLE SHEET"            )
  cat("\n ... ### ------------------------------------- \n\n" )
  
  # nf core style sample sheets must have following columns
  ## Strandneess and paired should preferrably be obtained from IEM, but as for now supplied as input variables
  ## group replicate fastq_1 fastq_2 strandness
  ## group: try to get group from sample sheet, if not available, then use Sample_Project
  ## replicate: use if presernt in sample sheet, else create one using Sample_Project
  ## For each Sample_Name processed, a suffiux is added as:
  ## _nS_R1_001.fastq.qz. The S-counter  will just be numeric order 1:nrow
  ## The L00X wil be differnent if multiple libraries, as for now 001 is the only one supported
  
  
  data_table.nfcore <- data_table.out
  
  ## create group column
  if(!("group" %in% colnames(data_table.nfcore))){
    cat("\n ... creating 'group' column from ¨Sample_Project¨")
    data_table.nfcore$group <- data_table.nfcore$Sample_Project
  }
  # generate replicate column
  if(!("replicate" %in% colnames(data_table.nfcore))){
    cat("\n ... creating 'replicate' column from 'group'")
    u <- sapply(unique(data_table.nfcore$group), function(x) length(which(data_table.nfcore$group==x)))
    rep.vec <- data_table.nfcore$group
    for(i in 1:length(u)){
      rep.vec[which(rep.vec==names(u)[i])] <- 1:u[i]
    }
    data_table.nfcore$replicate <- as.integer(rep.vec)
  }
  # generate batch column
  if(!("batch" %in% colnames(data_table.nfcore))){
    cat("\n ... creating 'batch' column from ¨Sample_Project¨")
    data_table.nfcore$batch <- data_table.nfcore$Sample_Project
  }
  
  
  fastq.sufifx <- "_001.fastq.qz"
  data_table.nfcore$fastq_1 <- paste(data_table.nfcore$Sample_Name, "_S",1:nrow(data_table.nfcore), "_R1", fastq.sufifx, sep="")
  data_table.nfcore$fastq_1_id_path <- paste(data_table.nfcore$Sample_ID, "/",data_table.nfcore$Sample_Name, "_S",1:nrow(data_table.nfcore), "_R1", fastq.sufifx, sep="")
  data_table.nfcore$fastq_2 <- paste(data_table.nfcore$Sample_Name, "_S",1:nrow(data_table.nfcore), "_R2", fastq.sufifx, sep="")
  data_table.nfcore$fastq_2_id_path <- paste(data_table.nfcore$Sample_ID, "/",data_table.nfcore$Sample_Name, "_S",1:nrow(data_table.nfcore), "_R2", fastq.sufifx, sep="")
  data_table.nfcore$bam <- paste0(data_table.nfcore$Sample_Name, opt$bam_suffix)
  data_table.nfcore$counts_matrix_id <- paste0(data_table.nfcore$Sample_Name, opt$count_mat_siuffix)
  
  ## re-order sheet
  data_table.nfcore <- dplyr::select(.data = data_table.nfcore, c(group, replicate, fastq_1, fastq_2, everything()))
  
  
  # If to save updated sample sheet to use in next run/itteratve process
  if(opt$output_nfcore){
    out.file.name <- paste(opt$sample_sheet, "_nfcore.csv", sep="")
    if(grepl("[.]csv", opt$sample_sheet)) out.file.name <- gsub(pattern = "[.]csv", replacement = "_nfcore.csv", x = opt$sample_sheet)
    cat("\n\n ... writing NF-core style sample sheet: ", out.file.name)
    write.table(x = data_table.nfcore, file = out.file.name, sep = ",", append = F, quote = F, row.names = F, col.names = T, na="", )
  }
  
  cat("\n ... ### --------------------------------### " )
  cat("\n ... ###   Wrote NF Core sample Sheet " )
  cat("\n ... ### ------------------------------- ### \n\n" )
}