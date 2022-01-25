# Creating a data dict for current UKBB application.
# Aims: Obtain table with all variable info (including sub-categories), re-categorising/ pre-processing informtion, naming guidance etc.

library(readr)
library(data.table)
library(openxlsx)

# Reading arguments
args <- commandArgs(trailingOnly = TRUE)
ukb_path <- as.character(args[1])

# Loading the data
ukb_main_columns <- colnames(fread(ukb_path, data.table = FALSE, nrows = 0)) # unique IDS from current UKBB basket

# Loading UK Biobank files
dict <- fread("../docs/Data_Dictionary_Showcase.csv", data.table = FALSE) # Data dictionary from UKBB: https://biobank.ctsu.ox.ac.uk/crystal/exinfo.cgi?src=accessing_data_guide
codings <- fread("../docs/Codings.csv", data.table = FALSE) # Codings for categorical outcomes: https://biobank.ctsu.ox.ac.uk/crystal/exinfo.cgi?src=accessing_data_guide

# Extracting field IDs
IDs <- unique(gsub("-.*", "", gsub("X", "", ukb_main_columns))) # unique IDS from application. Removing array & instance info.
matchedDict <- dict[dict$FieldID %in% IDs, ]
print("Field IDs that could not be matched:")
print(IDs[!IDs %in% intersect(IDs, unique(dict$FieldID))])

# Identifying columns to keep from UK Biobank showcase file
requiredColumns <- c("Field", "Participants", "Notes", "Units", "Path", "Link", "FieldID", "ValueType", "Coding") # Value type (integer, factor etc.) can be used to filter table for pre-processing sep differences between continious and cat. variables
matchedDict <- matchedDict[, requiredColumns] # Dictionary uncluding current application variables only

# Extra columns
matchedDict$Required <- NA # Is this variable required for analysis: no if NA, yes if anything else
matchedDict$RequiredCategorical <- NA # Is this variable required as categorical? #Options: 1 (Yes), 0 (No)
matchedDict$RequiredCatMethod <- NA # Categorising method: NA (nothing), 1 (Quartiles), 2 (Median), 3 (Mean), else put cut-off values in quotation marks, separated by a comma. Eg. age = "<50", "50-54", ">64"

matchedDict$GroupingName <- NA # Group assignment. e.g. socio-economic, psychosocial, environmental exp.
matchedDict$CodingName <- NA # shorter/ consice variable name for scripts - no white spaces etc.
matchedDict$FigureName <- NA # Names to put in figure - can include white spaces etc.
matchedDict$AddUnit <- NA # Indicates if unit should be added in the FigureName. If NA or FALSE, it is not added. If TRUE, the unit provided in "Unit" is added in brackets at the end of "FigureName". User-defined units can be specified as characters and will be used instead of the value provided in "Unit".

matchedDict$OG_categoriesNum <- NA # Sub-categories in numeric form - from coding showcase
matchedDict$OG_categoriesChar <- NA # Sub-categories in character form - from coding showcase

matchedDict$NewCatIDMatch <- NA # Numeric ID for re-categorisation of sub-categories. This ID pairs old sub-categories to new category assignment. See qualifications example below.
matchedDict$NewCatIDRef <- NA # ID of the category to use as reference
matchedDict$NewCatChar <- NA # New names for re-categorised variables, has to be in line with "NewCatIDMatch".
# eg. qualifications:  OG_categoriesNum = 1,2,3,4,5,6 ; OG_categoriesChar = GCSE, O-level, A-level, UG, PG, PhD ; NewCatIDMatch = 1,2,2,3,3,3 (n = 3) ; NewCatIDRef = 3 ; NewCatChar = Low, Medium, High (High is used as reference)
# to remove unknown/prefer not to answer for integer/continuous variables (e.g. number in household), OG_categoriesNum=-1,-3 ; OG_categoriesChar="Do not known", "Prefer not to answer" ; NewCatIDMatch=NA ; NewCatChar=NA (i.e. removed by default)
# to code unknown/prefer not to answer to be able to identify them in the data: OG_categoriesNum=-1,-3 ; OG_categoriesChar="Do not known", "Prefer not to answer" ; NewCatIDMatch=1,1 ; NewCatChar="exclude"
# the generated variable would become a character string

matchedDict$ArrayList <- NA
matchedDict$ArrayMethod <- NA # 1 = min, 2 = max, 3 = mean, 4 = median, 5 = first, 6 = last
# e.g. qualifications: to take the highest qualification across the multiple answers, we would use InstanceMethod = 1 (i.e. take the min NewCatIDMatch)

matchedDict$InstanceList <- NA
matchedDict$InstanceRequired <- NA # List of time points to keep. If NA, only baseline (0) is kept.
# e.g. To keep baseline and first repeat, use 0,1

# Defining ArrayList and InstanceList: lists of available arrays/instances
print("Defining available instances and arrays")
pb <- utils::txtProgressBar(style = 3)
for (i in 1:nrow(matchedDict)) {
  tmpname <- ukb_main_columns[grep(paste0("^", matchedDict[i, "FieldID"], "\\-"), ukb_main_columns)]
  tmpname <- gsub("-", ".", tmpname)
  matchedDict[i, "ArrayList"] <- paste(unique(sapply(strsplit(tmpname, split = "\\."), FUN = function(x) {
    x[3]
  })), collapse = ",")
  matchedDict[i, "InstanceList"] <- paste(unique(sapply(strsplit(tmpname, split = "\\."), FUN = function(x) {
    x[2]
  })), collapse = ",")
  utils::setTxtProgressBar(pb, i / nrow(matchedDict))
}
cat("\n")

# Adding information from original coding
print("Assigning coding IDs")
pb <- utils::txtProgressBar(style = 3)
for (i in 1:nrow(matchedDict)) {
  if (!is.na(matchedDict[i, "Coding"])) {
    codeID <- matchedDict[i, "Coding"]
    tmp <- codings[codings$Coding %in% codeID, ]
    matchedDict[i, "OG_categoriesNum"] <- paste0("\"", paste0(tmp$Value, collapse = "\",\""), "\"")
    matchedDict[i, "OG_categoriesChar"] <- paste0("\"", paste0(tmp$Meaning, collapse = "\",\""), "\"")
  }
  utils::setTxtProgressBar(pb, i / nrow(matchedDict))
}
cat("\n")

# Saving the prepared table
write_csv(matchedDict, "../docs/matchedDict.csv")
write.xlsx(matchedDict, "../docs/matchedDict.xlsx")
