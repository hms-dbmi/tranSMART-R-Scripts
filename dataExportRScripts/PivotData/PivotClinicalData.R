###########################################################################
 # tranSMART - translational medicine data mart
 #
 # Copyright 2008-2012 Janssen Research & Development, LLC.
 #
 # This product includes software developed at Janssen Research & Development, LLC.
 #
 # This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
 # as published by the Free Software  * Foundation, either version 3 of the License, or (at your option) any later version, along with the following terms:
 # 1.	You may convey a work based on this program in accordance with section 5, provided that you retain the above notices.
 # 2.	You may convey verbatim copies of this program code as you receive it, in any medium, provided that you retain the above notices.
 #
 # This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS    * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 #
 # You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.
 #
 #
 ##########################################################################


###########################################################################
#PivotClinicalData
#Parse the i2b2 output file and create input files for Cox/Survival Curve.
###########################################################################
Logger.debug <- function(msg) {
  print(paste(format(Sys.time(), "%Y-%m-%d %H:%M:%OS2"), msg))
}

PivotClinicalData.pivot <-function(inputFileName, snpDataExists, multipleStudies, study) {

library(data.table)
library(reshape2)

  Logger.debug(paste("Start reading input file", inputFileName))
  dataFile <- fread( input = inputFileName,
                        header = TRUE,
                        colClasses = "character")
  Logger.debug(paste("Finished reading input file", inputFileName))

  Logger.debug("Start sorting input data.")
  # Line below throws an error, easier just to copy the dataframe and work
  # on the copy.
  # subset <- dataFile[, c("PATIENT ID", "CONCEPT_PATH_FULL", "VALUE")]
  subset <- dataFile
  subset$pair <- paste( subset$`PATIENT ID`, subset$`CONCEPT_PATH_FULL`, sep="*")
  subset <- subset[!duplicated( subset$pair ), ]
  Logger.debug("Finished de-duping subset object.")

  Logger.debug("Start acast function run.")
  tt <- acast(subset, subset$`PATIENT ID`~subset$`CONCEPT_PATH_FULL`, value.var="VALUE", fill = " ")
  finalData <- as.data.frame( tt )
  finalData$PATIENT_ID <- rownames(finalData)
  finalData <- finalData[, c(ncol(finalData), 1:ncol(finalData)-1)]
  Logger.debug("Finished processing patient matrix.")

  filename <- "clinical_i2b2trans.csv"
  if (multipleStudies) filename <- paste(study, "_clinical_i2b2trans.csv")
  Logger.debug(paste("Start writing output file:", filename))
  write.csv(finalData, file = filename, row.names = FALSE)
  Logger.debug(paste("Finished writing output file."))

  Logger.debug("Done.")
}

PivotClinicalData.pivotDeprecated <-
function
(
input.dataFile, snpDataExists, multipleStudies, study
)
{
  print(snpDataExists)

	#Read the input file.
	dataFile <- data.frame(read.delim(input.dataFile))

	#Split the data by the CONCEPT_PATH.
	splitData <- split(dataFile,dataFile$CONCEPT.PATH)
  foo <- unique(dataFile[c("PATIENT.ID")])

  if (snpDataExists) {
    snpPEDFileData <- unique(subset(dataFile[c("PATIENT.ID", "SNP.PED.File")], SNP.PED.File != ""))
    colnames(snpPEDFileData) <- c("PATIENT.ID", "SNP.PED.File")
  }
	#Create a matrix with unique patient_nums.
	finalData <- matrix(unique(dataFile$PATIENT.ID));

  #Name the column.
	colnames(finalData) <- c("PATIENT.ID");

	#Get the unique list of concepts.
	conceptList <- unique(dataFile$CONCEPT.PATH)
	#For each of the passed in concepts, append the rows onto the end of our temp matrix.
	for(entry in conceptList)
	{
		#For each concept we merge the data against the patient num.
		finalData <- merge(finalData,unique(splitData[[entry]][c('PATIENT.ID','VALUE')]),by="PATIENT.ID",all.x=TRUE)
	}

  finalData <- merge(finalData, foo[c("PATIENT.ID")], by="PATIENT.ID", all.x=TRUE)

  if (snpDataExists) {
    finalData <- merge(finalData, snpPEDFileData, by="PATIENT.ID", all.x=TRUE)
    colnames(finalData) <- c("PATIENT ID",paste(conceptList), "SNP PED File")
  } else {
    colnames(finalData) <- c("PATIENT ID",paste(conceptList))
  }

  #We need MASS to dump the matrix to a file.
	require(MASS)

  # DI-811 Produce two similarly, yet differently formatted files
  # since the matrix option leaves padding in place. Remove padding
  # and wrap character cells in double-quotes, for further processing.
  filename <- "clinical_i2b2trans.txt"
  if (multipleStudies) filename <- paste(study, "_clinical_i2b2trans.txt")
  write.table(finalData, file = filename, row.names = FALSE, eol = "\n", sep = "\t")
  #file.remove(input.dataFile)
}
