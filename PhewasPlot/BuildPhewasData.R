#' PhewasData.build
#'
#' @export
PhewasData.build <- function(
  # input data -------
  input.dataFile = '~/SSC/transmartExtract/clinical.i2b2trans',
  output.file = 'phewas_data.json',
  output.metadata = 'metadata.json',
  # --------

  # groups variable -------
  group1Variable = "\\DBMI\\Autism\\SSC_wigler_mutations\\CHD8\\yes\\",
  group2Variable = "\\DBMI\\Autism\\SSC_wigler_mutations\\CHD8\\no\\",
  levelForCategories = 4,
  # --------

  # analysis parameters ------
  longitudinal = FALSE,
  maxUniqueValuesRatio = 0.9,

  scale =FALSE,
  logTransform = FALSE,

  threshold= 5,
  group_threshold = 2,

  discretization = FALSE,
  discMethod = 'cluster',
  ncategories = 3,

  ncore = 4,

  txtGLMFamily = 'binomial',
  numGLMFamily = 'gaussian'
  # ----------
  )
  {

  # source external functions ------
   #source('R/phewasPlot.R')
  # ------------

  # load external libraries --------
  library(reshape2)
  library(ggplot2)
  library(qqman)
  library(data.table)
  library(dplyr)
  library(arules)
  library(doSNOW)
  #library(rjson)
  library(ade4)
  library(jsonlite)
  # ----------
  #

        levelForCategories = as.integer(levelForCategories);
        maxUniqueValuesRatio = as.numeric(maxUniqueValuesRatio);
        ncategories = as.integer(ncategories);
        logTransform = as.logical(logTransform)
        scale = as.logical(scale)
        discretization = as.logical(discretization)

  print("-------------------")
  print("BuildPhewasData.r")

  print('Loading data')
  # load data -------
  data <- loadTranSMARTData(input.dataFile)
  # --------
  #
  #get patients groups from groupsvariables
  patientsExtrated <- getPatientGroupFromGroupVar(data, group1Variable, group2Variable)
  data <- patientsExtrated$data
  patientsGroups <- patientsExtrated$patientGroups
  rm(patientsExtrated)

  cl <- registerCluster(ncore,'SOCK')
  data <- pivotTransmartData(data,cl)
  stopCluster(cl)
  # levelForCategories <- levelForCategories - data$max_common_path_len
  # data <- data$dataCasted

  print('Pivot data')
  # pivot data --------
  preprocessed <- preprocessConcepts(data, longitudinal, maxUniqueValuesRatio)

  facts <- merge(patientsGroups,as.data.frame(preprocessed$facts), by='PATIENT_ID')

  if (sum(facts$group == 0) == 0 | sum(facts$group == 1) == 0) {
          stop('||FRIENDLY|| There is no data for one of the groups' )
  }

  concepts <- preprocessed$concepts
  rm(preprocessed)
  # ----------

  print('Transform numerical concepts')
  # transform num concepts -------------
  numConcepts <- concepts$CONCEPT_PATH[concepts$type == 'num']

  if (scale) {
    scaled <- scaleConcepts(facts,numConcepts)

    facts <- scaled$facts
    scaleParams <- scaled$scaleParams

    rm(scaled)
  }

  if (logTransform) {
    facts <- logTransformConcepts(facts,numConcepts)
  }

  if (discretization) {
    discretized <- discetizeConcepts(facts, concepts, numConcepts, discMethod, ncategories)
    facts <- discretized$facts
    concepts <- discretized$concepts
    rm(discretized)
  }
  # ----------
  #
  if (nrow(concepts) < ncore) {
          ncore <- 1
  }

  cl <- registerCluster(ncore,'SOCK')

  print('analysis loop')
  results <- analyzeConcepts(facts = facts,
                   concepts = concepts,
                   group_var = 'group',
                   covars=c(''),
                   threshold= threshold,
                   group_threshold = group_threshold,
                   scale = scale,
                   scaleParams = if (scale) {scaleParams} else {""} ,
                   logTransform = logTransform,
                   cluster = cl,
                   txtGLMfamily= txtGLMfamily,
                   numGLMfamily = numGLMfamily,
                   levelcat = levelForCategories)

  stopCluster(cl)
  print('write results table')

  metadata <- extractMetadata(group1Variable,
                  group2Variable,
                  levelForCategories,
                  longitudinal,
                  maxUniqueValuesRatio,
                  scale,
                  logTransform,
                  discretization,
                  discMethod,
                  ncategories,
                  ncore,
                  txtGLMFamily,
                  numGLMFamily)

#   write.table(results,file=output.file,col.names=T, row.names=F,append=F,sep='\t')

  #results$CONCEPT_PATH <- gsub('[\\]','>',results$CONCEPT_PATH)
  #results$CONCEPT_PATH_FULL <- gsub('[\\]','>',results$CONCEPT_PATH_FULL)

  cat(toJSON(results), file=output.file)
  cat(toJSON(metadata), file=output.metadata)
  print("-------------------")
  return(list(results= results, metadata = metadata))
}

#' loadData
#'
#' @export
loadData <- function(input.dataFile, group1Variable, group2Variable) {
  # load data -------
  if (grepl('\\.gz',input.dataFile)) {
          gunzip(paste0(input.dataFile))
          input.dataFile <- gsub('.gz','',input.dataFile)
  }
  data <- fread(input.dataFile, stringsAsFactors = F)
  # ---------

  # remove spaces in column names ------
  setnames(data,names(data),gsub(' ','_',names(data)))
  # --------

  # filter facts for patients with known groups ---------
  group1 <- grepl(group1Variable,data$CONCEPT_PATH_FULL,fixed=T)
  group2 <- grepl(group2Variable,data$CONCEPT_PATH_FULL,fixed=T)
  group1Patients <- data$PATIENT_ID[group1]
  group2Patients <-  data$PATIENT_ID[group2]
  data <- data[!group1 & !group2,]
  data <- data[data$PATIENT_ID %in% c(group1Patients,group2Patients),]
  # ---------

  # generate group variable ---------
  # group = 1 for patients in group 1 and 0 for patients in group 2 --------
  data$group <- 0
  data$group[data$PATIENT_ID %in% group1Patients] <- 1
  # -----------

  if (sum(data$group) == 0) {
        stop('There is no patient in group 1')
  } else if (sum(data$group) == nrow(data)) {
        stop('There is no patient in group 2 ')
  }

  # copy concept_path_full to concept_path to solve data extraction issues -------
  data$CONCEPT_PATH <- data$CONCEPT_PATH_FULL
  # ---------

  return(data)
}

#' loadTranSMARTData
#'
#' @export
loadTranSMARTData <- function(input.dataFile) {
        require(data.table)
        # load data -------
        if (grepl('\\.gz',input.dataFile)) {
                gunzip(paste0(input.dataFile))
                input.dataFile <- gsub('.gz','',input.dataFile)
        }
        data <- fread(input.dataFile, stringsAsFactors = F)
        # ---------

        # remove spaces in column names ------
        setnames(data,names(data),gsub(' ','_',names(data)))
        # --------

        return(data)
}

#' getPatientGroupFromGroupVar
#' @description get the list of patients and their analysis group from the extracted tranSMART data
#' @param data : data.frame with data extracted from tranSMART
#' @param group1Variable : variable defining the case group
#' @param group2Variable : variable defining the control group
#' @return list containing data without the grouping variables and patientGroups
#' @export
getPatientGroupFromGroupVar <- function(data, group1Variable, group2Variable) {

        # assign patients to groups from data ------
        group1 <- data.frame(PATIENT_ID = data$PATIENT_ID[data$CONCEPT_PATH_FULL == group1Variable], group = 1, stringsAsFactors = F)
        group2 <- data.frame(PATIENT_ID = data$PATIENT_ID[data$CONCEPT_PATH_FULL == group2Variable], group = 0, stringsAsFactors = F)

        # if one of the groups is empty -> error -----
        if (nrow(group1) == 0) {
                stop('||FRIENDLY||There is no patients in group 1')
        } else if (nrow(group2) == 0) {
                stop('||FRIENDLY||There is no patients in group 2')
        }
        # ---------

        # deduplicate patient list -------
        patientGroups <- rbind(group1,group2)
        patientGroups <- unique(patientGroups)
        # --------

        # remove grouping variable from data -------
        data <- data[(!data$CONCEPT_PATH_FULL %in% c(group1Variable,group2Variable)),]
        # --------

        return(list(data = data, patientGroups = patientGroups) )

}

#' registerCluster
#' @description registers a cluster using doSNOW
#' @param ncore: number of cores
#' @param type: cluster type
#' @return cluster object
#' @export
registerCluster <- function(ncore,type) {
        require(foreach)
        require(doSNOW)
        require(data.table)
        cl <- makeCluster(ncore, type)
        registerDoSNOW(cl)
        return(cl)
}

#' LCS
#' @description computes the longuest common sequence from 2 string vectors.
#' @param s1 : string vector
#' @param s2 : string vector
#' @return longuest common vector
#' @export
LCS <- function(s1,s2) {
        l <- min(length(s1),length(s2))
        for (i in 1:l) {
                if (s1[i] != s2[i]) {
                        return (s1[1:(i-1)])
                }
        }
        return (s1[1:i])
}

#' getMaxCommonPathFromConceptPath
#' @description computes the longuest common sequence from a list of concept_paths
#' @param concept_paths : vector of concept_paths (point delimited)
#' @return list of max common path and length of this path
#' @export
getMaxCommonPathFromConceptPath <- function(concept_paths) {
        concept_paths <- unique(concept_paths)
        concept_paths <- strsplit(concept_paths, '\\.')

        max_common_path <- concept_paths[[1]]

        if (length(concept_paths) > 1) {
                for (i in 2:length(concept_paths)) {
                        max_common_path <- LCS(max_common_path, concept_paths[[i]])
                }
        }

        len <- length(max_common_path)

        max_common_path <- paste(max_common_path, collapse = '.')
        max_common_path <- paste0(max_common_path,'.')
        return(list(max_common_path = max_common_path, len = len))
}

#' pivotTransmartData
#' @description transforms tranSMART data from long to wide format. i.e. 1 obvervation per patient and one variable per concept
#' @param data : extracted data from tranSMART
#' @param cl: cluster
#' @return data.frame of pivoted data
#' @export
pivotTransmartData <- function(data, cl) {

        data$characterValue<- sapply(data$VALUE, function(x) {
                is.na(as.numeric(x))
        })

        data$CONCEPT_PATH_FULL <- gsub('^[\\]','',data$CONCEPT_PATH_FULL)
        data$CONCEPT_PATH_FULL <- gsub('[\\]','\\.',data$CONCEPT_PATH_FULL)
        max_common_path <- getMaxCommonPathFromConceptPath(data$CONCEPT_PATH_FULL)
        max_common_path_len <- max_common_path$len
        max_common_path <- max_common_path$max_common_path



        data$CONCEPT_PATH_FULL <- gsub(max_common_path,'',data$CONCEPT_PATH_FULL)

        patients <- data$PATIENT_ID[!duplicated(data$PATIENT_ID)]
        patients <- data.frame(PATIENT_ID=patients[order(patients)],stringsAsFactors = F)

        data$CONCEPT_PATH[data$characterValue] <- gsub('[.][^.]+[.]$','',data$CONCEPT_PATH_FULL[data$characterValue])
        #data$CONCEPT_PATH[!data$characterValue] <- paste0(gsub('[\\][^\\]+[\\]$','',data$CONCEPT_PATH_FULL[!data$characterValue]),'.num')
        data$CONCEPT_PATH[!data$characterValue] <- paste0(data$CONCEPT_PATH_FULL[!data$characterValue],'num')
        data$CONCEPT_PATH <- paste0(max_common_path,data$CONCEPT_PATH)

        testdata <- split(data,data$CONCEPT_PATH)

        recombine <- function(concept,patients) {
                casted <- dcast.data.table(concept, PATIENT_ID ~ CONCEPT_PATH, value.var = 'VALUE')
                casted <- merge(patients,casted, by='PATIENT_ID', all.x=T)
                return(casted)
        }

        dataCasted <- foreach(x = 1:length(testdata), .export =c('patients'), .packages = 'data.table', .combine= cbind) %do%
                recombine(testdata[[x]],patients)

        IDCols <- names(dataCasted) == 'PATIENT_ID'
        IDCols[1] <- FALSE

        dataCasted <- dataCasted[, !IDCols]

        return(dataCasted)
}


#' preprocessConcepts
#'
#' @export
preprocessConcepts <- function(data, longitudinal, maxUniqueValuesRatio) {
        # manage duplicated facts ---------
        if (longitudinal) {
                #TODO: manage longitudinal data with aggregation functions: oldest, most recent, mean, sum, max, min, others...
        } else {
        #        setkeyv(data, c('PATIENT_ID','CONCEPT_PATH'))
                data <- unique(data)
        }
        # --------

        names(data) <- gsub('[\\]', '.',names(data))
        names(data) <- gsub('^\\.', '',names(data))

        preproFacts <- subset(data, select = 'PATIENT_ID')
        preproConcepts <- data.frame(data.frame(CONCEPT_PATH = character(), levels =  numeric(), NAs = numeric(),
                                                count = numeric(), NonZeroCount = numeric(), uniqueValueRatio = numeric(),
                                                type = character(), stringsAsFactors = F)
        )

        getConceptsStats <- function(data) {
                levels <- sapply(2:ncol(data), function(x) length(levels(as.factor(data[,x]))))
                nas <- sapply(2:ncol(data), function(x) sum(is.na(data[,x])))
                count <- sapply(2:ncol(data), function(x) sum(!is.na(data[,x])))
                NonZeroCount <- sapply(2:ncol(data), function(x) sum(!is.na(data[,x]) & data[,x] != 0))
                concepts <- data.frame(CONCEPT_PATH = names(data)[-1], levels =  levels, NAs = nas, count = count, NonZeroCount, stringsAsFactors = F)
                concepts$uniqueValueRatio <- concepts$levels / concepts$count
                return(concepts)
        }

        concepts <- getConceptsStats(data)
        concepts$filterOut <- concepts$levels == 1 | concepts$uniqueValueRatio >= maxUniqueValuesRatio

        data <- data[,names(data) %in% concepts$CONCEPT_PATH[!concepts$filterOut] | names(data) == 'PATIENT_ID' ,]
        concepts <- concepts[!concepts$filterOut,]
        concepts$filterOut <- c()

        if (nrow(concepts) == 0) { stop('||FRIENDLY|| No concept could be analyzed according to the selected analysis parameters')}

        concepts$character <-  sapply(2:ncol(data), function(x)  sum(is.na(as.numeric(data[!is.na(data[,x]),x]))) > 0)

        dataTxt <- subset(data, select= names(data) %in% concepts$CONCEPT_PATH[concepts$character & concepts$levels > 2])
        if (ncol(dataTxt) > 0) {
                dataTxt$BLANK <- 0
                dataTxtBin <- acm.disjonctif(dataTxt)
                dataTxtBin <- lapply(1:ncol(dataTxt), function(x) {
                        varX <- grep(names(dataTxt)[x], names(dataTxtBin))
                        dataTxtBin[is.na(dataTxt[,x]),varX] <- NA
                        dataTxtBin[,varX]
                })
                dataTxtBin <- data.frame(dataTxtBin)
                dataTxtBin$BLANK.0 <- c()
                dataTxtBin <- cbind(PATIENT_ID = data$PATIENT_ID, dataTxtBin, stringsAsFactors = F)
                conceptsTxtBin <- getConceptsStats(dataTxtBin)
                conceptsTxtBin$type = 'bin'

                preproFacts <- merge(preproFacts, dataTxtBin, all.x=T, by = 'PATIENT_ID')
                preproConcepts <- rbind(preproConcepts, conceptsTxtBin)
        }


        dataBin <- subset(data, select=(names(data) %in% concepts$CONCEPT_PATH[concepts$levels == 2]))
        if(ncol(dataBin) > 0) {
                dataBin01 <- sapply(1:ncol(dataBin), function(x) {
                        dd <- dataBin[,x]
                        levels <- data.frame(table(as.factor(dd)))
                        minLevel <- levels$Var1[levels$Freq == min(levels$Freq)][1]
                        dd[dd != minLevel & !is.na(dd)] <- 0
                        dd[dd == minLevel & !is.na(dd)] <- 1
                        dd <- as.numeric(dd)
                        dd <- list(dd)
                        names(dd) <- paste0(names(dataBin)[x],'.',minLevel)
                        return (dd)
                })
                dataBin <- data.frame(dataBin01)
                dataBin <- cbind(PATIENT_ID = data$PATIENT_ID, dataBin, stringsAsFactors = F)
                conceptsBin <- getConceptsStats(dataBin)
                conceptsBin$type = 'bin'

                preproFacts <- merge(preproFacts, dataBin, all.x=T, by = 'PATIENT_ID')
                preproConcepts <- rbind(preproConcepts, conceptsBin)
        }

        dataNum <- subset(data, select=(names(data) %in% concepts$CONCEPT_PATH[!concepts$character & concepts$levels > 2] ))
        if (ncol(dataNum) > 0) {
                dataNum <-sapply(1:ncol(dataNum), function(x) {
                        dd <- dataNum[,x]
                        dd <- as.numeric(dd)
                        dd <- list(dd)
                        names(dd) <- names(dataNum)[x]
                        return(dd)
                })
                dataNum <- data.frame(dataNum)
                dataNum <- cbind(PATIENT_ID = data$PATIENT_ID, dataNum, stringsAsFactors = F)
                conceptsNum <- concepts[!concepts$character & concepts$levels > 2,]
                conceptsNum$character <- c()
                conceptsNum$type = 'num'

                preproFacts <- merge(preproFacts, dataNum, all.x=T, by = 'PATIENT_ID')
                preproConcepts <- rbind(preproConcepts, conceptsNum)
        }


        #data <- merge(dataBin, dataNum, by = 'PATIENT_ID', all.x = T, all.y = T)
        #data <- merge(data, dataTxtBin,  by = 'PATIENT_ID', all.x = T, all.y = T)

        #concepts <- rbind(conceptsBin, conceptsNum, conceptsTxtBinX)

        preproConcepts$filterOut <- preproConcepts$levels == 1 | preproConcepts$uniqueValueRatio >= maxUniqueValuesRatio
        data <- data[,names(data) %in% preproConcepts$CONCEPT_PATH[!preproConcepts$filterOut] | names(data) == 'PATIENT_ID' ,]
        preproConcepts <- preproConcepts[!preproConcepts$filterOut,]
        preproConcepts$filterOut <- c()

        return(list(facts = preproFacts, concepts = preproConcepts))

}

preprocessConcepts_bak <- function(data, longitudinal, maxUniqueValuesRatio) {
        # manage duplicated facts ---------
        if (longitudinal) {
                #TODO: manage longitudinal data with aggregation functions: oldest, most recent, mean, sum, max, min, others...
        } else {
                setkeyv(data, c('PATIENT_ID','CONCEPT_PATH'))
                data <- unique(data)
        }
        # --------

        # count number of occurences of concepts ---------
        concepts <- subset(data,select=c('CONCEPT_CODE','CONCEPT_PATH','CONCEPT_PATH_FULL','VALUE'))
        setkeyv(concepts, c('CONCEPT_CODE','CONCEPT_PATH_FULL','VALUE'))

        concepts$ValueCount <- 1
        concepts <- aggregate(ValueCount ~ CONCEPT_CODE + CONCEPT_PATH + CONCEPT_PATH_FULL + VALUE, data = concepts, FUN = sum)

        concepts <- unique(concepts)
        # -----------

        # find concept type -------------
        concepts$type <- 'text'
        concepts$count <- 1

        # find numeric concepts ----------
        # if the end of CONCEPT_PATH == end of CONCEPT_PATH_FULL modulo '\\' then it's a numeric concept.

        m <- regexpr('[^\\]+$',concepts$CONCEPT_PATH)
        pathEnd <- regmatches(concepts$CONCEPT_PATH, m)

        m <- regexpr('[^\\]+[\\]$',concepts$CONCEPT_PATH_FULL)
        pathFullEnd <- regmatches(concepts$CONCEPT_PATH_FULL, m)
        pathFullEnd <- gsub('[\\]$','',pathFullEnd)

        concepts$End <- pathFullEnd


        concepts$type[!is.na(as.numeric(concepts$VALUE))] <- 'num'

        conceptsLevels <- aggregate(count ~ CONCEPT_PATH  + type, data = concepts, FUN = sum )

        # ----------

        # remove concepts with only 1 level --------
        concepts <- concepts[concepts$CONCEPT_PATH %in% conceptsLevels$CONCEPT_PATH[conceptsLevels$count !=1],]
        conceptsLevels <- conceptsLevels[conceptsLevels$count > 1,]
        # ----------

        conceptGroup = 1
        concepts$conceptGroup = 0
        toRemove <- c()

        for (i in 1:nrow(conceptsLevels)) {
                if (conceptsLevels$type[i] == 'text') {
                        c <- concepts[concepts$CONCEPT_PATH == conceptsLevels$CONCEPT_PATH[i],]
                        concepts <- concepts[concepts$CONCEPT_PATH != conceptsLevels$CONCEPT_PATH[i],]

                        if (conceptsLevels$count[i] == 2) {

                                c <- c[order(c$ValueCount),]
                                toRemove <- c(toRemove,c$CONCEPT_PATH_FULL[2])
                                c$CONCEPT_PATH = paste0(c$CONCEPT_PATH[1],'\\',c$VALUE[1])
                                c$VALUE = c(1,0)

                        } else {
                                c$conceptGroup = conceptGroup
                                conceptGroup = conceptGroup + 1

                                c$CONCEPT_PATH <- paste0(c$CONCEPT_PATH,'\\',c$VALUE)
                                c$VALUE = 1
                        }
                        if (conceptsLevels$count[i] /sum(c$ValueCount) <= maxUniqueValuesRatio) {
                                concepts <- rbind(concepts,c)
                        }
                }
        }

        concepts <- concepts[!duplicated(concepts$CONCEPT_CODE),]

        data <- data[data$CONCEPT_CODE %in% concepts$CONCEPT_CODE,]

        dataTxt <- merge(subset(data,select=c('CONCEPT_CODE',  'PATIENT_ID',	'SUBSET','group')),
                         subset(concepts[concepts$type == 'text',],select=c('CONCEPT_CODE','CONCEPT_PATH','CONCEPT_PATH_FULL','VALUE','conceptGroup')), by = 'CONCEPT_CODE')

        dataNum <- merge(subset(data,select=c('CONCEPT_CODE',  'PATIENT_ID',  'SUBSET','group', 'VALUE')),
                         subset(concepts[concepts$type == 'num',],select=c('CONCEPT_CODE','CONCEPT_PATH','CONCEPT_PATH_FULL','conceptGroup')), by = 'CONCEPT_CODE')


        # change the concept path depending on the type and number of levels -------
        # type == num: path doesn't change
        # type == text and count == 2 : concept_path = concept_path_full of the minority class and value = 1 for minority class
        # type == text and count > 2: concept_path = concept_path_full


        pivoted <- data[!duplicated(data$PATIENT_ID), c('PATIENT_ID','group'), with = F]

        if (nrow(dataNum) > 0) {
                dataNum$VALUE = as.numeric(dataNum$VALUE)
                dataNum <- dcast.data.table(data = dataNum,
                                            formula = PATIENT_ID  ~ CONCEPT_PATH,
                                            value.var = 'VALUE',
                                            fun.aggregate = sum,
                                            fill = NA)

                pivoted <- merge(pivoted, dataNum, by = 'PATIENT_ID', all.x = T)
        }


        conceptGroups <- levels(as.factor(concepts$conceptGroup))
        #conceptGroups <- conceptGroups[conceptGroups != 0]


        if (nrow(dataTxt) > 0) {
                dataTxt$VALUE <- as.numeric(dataTxt$VALUE)

                for (i in conceptGroups) {

                        temp <- dcast.data.table(dataTxt[dataTxt$conceptGroup == i,],
                                                 formula = PATIENT_ID  ~ CONCEPT_PATH,
                                                 value.var = 'VALUE',
                                                 fun.aggregate = sum,
                                                 fill = 0)
                        pivoted <- merge(pivoted, temp, by = 'PATIENT_ID', all.x = T)
                }
        }

        concepts <- concepts[concepts$CONCEPT_PATH %in% names(pivoted),c('CONCEPT_CODE','CONCEPT_PATH','CONCEPT_PATH_FULL','type','End','conceptGroup')]
        concepts <- concepts[!(concepts$CONCEPT_PATH_FULL %in% toRemove),]

        concepts$CONCEPT_CODE <- paste0('c',concepts$CONCEPT_CODE)

        setnames(pivoted, concepts$CONCEPT_PATH, concepts$CONCEPT_CODE)

        return(list(facts = pivoted, concepts = concepts))

}


#' scaleConcepts
#'
#' @export
scaleConcepts <- function(facts, numConcepts) {

  sds <- unlist(lapply(numConcepts, function(x) {
    sd(facts[,x],na.rm=T)
  }))

  scaleParams <- data.frame(concept= numConcepts, sd = sds, stringsAsFactors = F)

  for (col in numConcepts) {
    facts[,col] <- (facts[,col])/sd(facts[,col],na.rm = T)
  }
  return(list(facts = facts, scaleParams = scaleParams))
}

#' logTransformConcepts
#'
#' @export
logTransformConcepts <- function(facts, numConcepts) {
  for (col in numConcepts) {
    facts[,col] <- log(facts[,col] + 1)
  }
  return(facts)
}

#' discretizeConcepts
#'
#' @export
discetizeConcepts <- function(facts, concepts,numConcepts, discMethod, ncategories ) {
  n <- length(numConcepts)
  for (i in 1:n) {
    colName <- names(facts)[names(facts) == numConcepts[i]]
    cat(paste0('[',i,'/',n,'] ',colName,'...'))
    if (sum(!is.na(facts[,numConcepts[i]])) >= ncategories) {
      col <- discretize(facts[,numConcepts[i]], method='interval', categories = ncategories)
    } else {
      col <- discretize(facts[,numConcepts[i]], method=discMethod, categories = ncategories)
    }
    newCols <- as.data.frame(lapply(1:ncategories, function(x) {as.numeric(col == levels(col)[x]) }))
    names(newCols) <- paste0(colName,'_',levels(col))
    con <- concepts[concepts$CONCEPT_PATH == colName,]
  newCons<- data.frame(CONCEPT_CODE = paste0(con$CONCEPT_PATH,'_',1:ncategories),
                       CONCEPT_PATH = names(newCols),
                       CONCEPT_PATH_FULL = paste0(con$CONCEPT_PATH_FULL,'\\',levels(col)),
                       End = paste0(con$End,'_',levels(col)),
                       type = 'text',
                       conceptGroup = 0)
    facts[,colName] <- c()
    concepts <- concepts[concepts$CONCEPT_PATH != colName,]
    facts <- cbind(facts,newCols)
    concepts <- rbind(concepts, newCons)
    cat('Done\n')
  }
  return(list(facts = facts, concepts= concepts))
}

#' analyzeConcepts
#'
#' @export
analyzeConcepts <- function(facts,concepts ,group_var = '', covars='',threshold= 5, group_threshold = 0, logTransform = logTransform,scale = scale, scaleParams = '', cluster = cl, txtGLMfamily= 'binomial', numGLMfamily = 'gaussian',
                            levelcat = levelForCategories) {

  # set groups and detect if the number of groups != 2 ------
  group_level <- levels(as.factor(facts[,group_var]))
  print(paste0('case group set to ',group_level[2]))
  if(length(group_level) != 2) (stop('Number of groups for group_var must be equal to 2'))
  # ---------

  # calculate NAs, number of cases and number of controls available for each concept
  concepts$NAs <- unlist(parLapply(cluster,concepts$CONCEPT_PATH, function(x) (sum(is.na(facts[,x])))))
  concepts$ncases <- unlist(parLapply(cluster,concepts$CONCEPT_PATH, function(x) (sum(!(is.na(facts[facts[,group_var]==group_level[2],x])) & facts[facts[,group_var]==group_level[2],x] != 0))))
  concepts$ncontrols <- unlist(parLapply(cluster,concepts$CONCEPT_PATH, function(x) (sum(!(is.na(facts[facts[,group_var]==group_level[1],x])) & facts[facts[,group_var]==group_level[1],x] != 0))))
  concepts$mean_cases <- unlist(parLapply(cluster,concepts$CONCEPT_PATH, function(x) (mean(facts[facts[,group_var]==group_level[2],x],na.rm=T))))
  concepts$mean_controls <- unlist(parLapply(cluster,concepts$CONCEPT_PATH, function(x) (mean(facts[facts[,group_var]==group_level[1],x],na.rm=T))))
  concepts$count <- concepts$ncases + concepts$ncontrols
  concepts$method <- ifelse(concepts$type == 'bin', 'glm','lm')
  # ---------

  # concepts without enough observations ------
  concepts <- concepts[(concepts$ncases + concepts$ncontrols) > threshold &
                         concepts$ncases > group_threshold &
                         concepts$ncontrols > group_threshold,]
  # ---------

  # Analysis loop --------
  pval <- loopGlm(data= facts,
                   x = group_var,
                   concept_list = concepts$CONCEPT_PATH,
                   type= concepts$type,
                   covar = covars ,cluster = cluster  )
  # --------

  if (!is.null(pval)) {
          concepts <- cbind(concepts, pval)

          if (logTransform) {
                  for (concept in concepts$CONCEPT_PATH[concepts$type == 'num']) {
                          concepts$mean_cases[concepts$CONCEPT_PATH == concept] <- exp(concepts$mean_cases[concepts$CONCEPT_PATH == concept]) - 1
                          concepts$mean_controls[concepts$CONCEPT_PATH == concept] <- exp(concepts$mean_controls[concepts$CONCEPT_PATH == concept]) - 1
                  }
          }

          if (scale) {
                  for(i in 1:nrow(scaleParams)) {
                          concepts$mean_cases[concepts$CONCEPT_PATH == scaleParams$concept[i]] <- (concepts$mean_cases[concepts$CONCEPT_PATH == scaleParams$concept[i]] * scaleParams$sd[i])
                          concepts$mean_controls[concepts$CONCEPT_PATH == scaleParams$concept[i]] <- (concepts$mean_controls[concepts$CONCEPT_PATH == scaleParams$concept[i]] * scaleParams$sd[i])
                  }
          }

          concepts$p_bonf <- p.adjust( concepts$pvalue, "bonferroni")
          concepts$p_fdr <- p.adjust(concepts$pvalue, "BH")


          concepts$category <- unlist(lapply(concepts$CONCEPT_PATH,function(x) {
                  res <- unlist(strsplit(x, '\\.'))
                  return(res[[levelcat]])
          }))

          concepts <- concepts[order(concepts$pvalue),]
          concepts <- concepts[!is.na(concepts$pvalue),]
          return(concepts)

  } else {
          stop('||FRIENDLY|| No concept could be analyzed according to the selected analysis parameters')
  }

}

#' loopGlm
#'
#' @export
loopGlm <- function (data, x, concept_list,type, covar = '', cluster = cluster, txtGLMfamily= 'binomial', numGLMfamily = 'gaussian' )
{

  if (length(concept_list) > 0) {
          pval <- foreach(i=1:length(concept_list), .combine='rbind', .export=c('runGlm')) %dopar%
          {
                  if (type[i] == 'text') {
                          runGlm(data=data, x=concept_list[i], y=x, covar=covar, GLMfamily = txtGLMfamily)
                  } else {
                          runGlm(data=data, x=concept_list[i], y=x, covar=covar, GLMfamily = numGLMfamily)
                  }
          }
  } else {
          pval <- NULL
  }

return(pval)
}

#' runGlm
#'
#' @export
runGlm <- function(data, x, y ,covar = '', GLMfamily = 'binomial' )
{
  if (covar[1] != '') {paste('as.numeric(',y,')+', paste(covar,collapse='+')) -> form.right} else {y -> form.right}

  if (GLMfamily == 'binomial') {
    family = binomial
  } else if (GLMfamily == 'gaussian') {
    family = gaussian
  }

  form <- as.formula(paste(x, '~', form.right))
  fit <- glm(form, data=data, family=family)
  coef <- exp(coef(fit)[2])
  CI <- exp(confint(fit)[2,])
#   if (GLMfamily == 'binomial') {
#     coef <- exp(coef)
#     CI <- exp(CI)
#   }
  coefs <- data.frame(coef = coef, CI_inf = CI[1], CI_sup = CI[2], pvalue = (summary(fit)$coef)[2,4])
  return(coefs)
}

#' toJSONarray
#'
#' @export
toJSONarray <- function(dtf){
  clnms <- colnames(dtf)
  name.value <- function(i){
    quote <- '';
    if(!class(dtf[, i]) %in% c('numeric', 'integer')){
      quote <- '"';
    }
    paste('"', i, '" : ', quote, dtf[,i], quote, sep='')
  }
  objs <- apply(sapply(clnms, name.value), 1, function(x){paste(x, collapse=', ')})
  objs <- paste('{', objs, '}')
  res <- paste('[', paste(objs, collapse=', '), ']')
  return(res)
}

#' extractMetadata
#'
#' @export
extractMetadata <- function(group1Variable,
                            group2Variable,
                            levelForCategories,
                            longitudinal,
                            maxUniqueValuesRatio,
                            scale,
                            logTransform,
                            discretization,
                            discMethod,
                            ncategories,
                            ncore,
                            txtGLMFamily,
                            numGLMFamily) {

  metadata <- data.frame(group1Variable = group1Variable,
                   group2Variable = group2Variable,
                   levelForCategories = levelForCategories,
                   longitudinal = longitudinal,
                   maxUniqueValuesRatio = maxUniqueValuesRatio,
                   scale =scale,
                   logTransform = logTransform,
                   discretization = discretization,
                   discMethod = discMethod,
                   ncategories = ncategories,
                   ncore = ncore,
                   txtGLMFamily = txtGLMFamily,
                   numGLMFamily = numGLMFamily)

  metadata$group1Variable <- gsub('[\\]','>',metadata$group1Variable)
  metadata$group2Variable <- gsub('[\\]','>',metadata$group2Variable)
  return(metadata)
}

#' return 1
#'
#' @export
return1 <- function() {
        return(1)
}

