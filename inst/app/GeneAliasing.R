
# Perform the actual query against STRINGdb
.QuerySTRINGdb <- function(inputIds, speciesId, score_threshold = 0, stringDBVersion = "10"){
	print('Querying STRINGdb')
  string_db <- STRINGdb::STRINGdb$new(version = stringDBVersion,
                            species = speciesId,
                            score_threshold = score_threshold,
                            input_directory = "")

  ## map inputIds to stringID mapIds
  inputIds.map <- string_db$map(my_data_frame = data.frame(InputTerm = as.character(inputIds)),
                                my_data_frame_id_col_names = "InputTerm",
                                takeFirst = T,
                                removeUnmappedRows = FALSE)
	## get all species aliases
  stringdb.alias <- string_db$get_aliases()
  stringdb.alias <- stringdb.alias %>%
    group_by(STRING_id) %>%
    summarize(STRING.aliases = toString(sort(unique(alias))))

  stringdb.alias <- merge(inputIds.map, stringdb.alias, by = c("STRING_id"), all.x = F)

	print(paste0('Found ', sum(!is.na(stringdb.alias$STRING_id)), ' of ', length(inputIds)))

  return(stringdb.alias)
}



# Basic argument checking
.CheckGeneInputs <- function(ensemblIds, geneSymbols){
  if (.IsEmpty(ensemblIds) && .IsEmpty(geneSymbols)) {
    stop('Must provide either ensemblIds or geneSymbols')
  }
  else if (!.IsEmpty(ensemblIds) & !.IsEmpty(geneSymbols)) {
    if (length(ensemblIds) != length(geneSymbols)) {
      stop('EnsemblIds and geneSymbol must be of equal length.')
    }
  }
}


# Utility function to test if the input is NA or NULL
.IsEmpty <- function(x) {
	return(all(is.na(x)) || all(is.null(x)))
}


#' @title TranslateToEnsembl
#' @param ensemblIds A vector of ensembl IDs, passed to the biomaRt::getBM() to query against ensembl_gene_id
#' @param geneSymbols A vector of gene symbols, passed to the biomaRt::getBM() to query against hgnc_symbol
#' @param dataset Passed directly to biomaRt::useEnsembl
#' @param ensemblVersion Passed directly to biomaRt::useEnsembl
#' @param ensemblMirror Passed directly to biomaRt::useEnsembl
#' @param biomart Passed directly to biomaRt::useEnsembl
#' @importFrom biomaRt useEnsembl getBM
#' @importFrom dplyr %>% group_by summarize mutate_all coalesce
#' @importFrom BiocGenerics order
#' @importFrom stringr str_match
#' @importFrom stats setNames
TranslateToEnsembl <- function(ensemblIds = NULL, geneSymbols = NULL, dataset = "mmulatta_gene_ensembl", ensemblVersion = NULL, ensemblMirror = "uswest", biomart = "ensembl"){
  .CheckGeneInputs(ensemblIds = ensemblIds, geneSymbols = geneSymbols)
  if (is.null(ensemblIds)) {
    ensemblIds <- NA
  }

  if (is.null(geneSymbols)) {
    geneSymbols <- NA
  }

	combinedEnsembl <- data.frame(EnsemblId = ensemblIds, GeneSymbol = geneSymbols, stringsAsFactors = FALSE)
	combinedEnsembl$Order <- 1:nrow(combinedEnsembl)

	ensemblById <- NA
	if (!.IsEmpty(ensemblIds)) {
		ensemblById <- .QueryEnsembl(inputIds = ensemblIds,
                             queryField = "ensembl_gene_id",
                             dataset = dataset,
                             ensemblVersion = ensemblVersion,
														 ensemblMirror = ensemblMirror,
                             biomart = biomart)

		ensemblById$EnsemblId <- ensemblById$ensembl_gene_id
		ensemblById <- merge(combinedEnsembl, ensemblById, by = 'EnsemblId', all.x = T)
		ensemblById <- dplyr::arrange(ensemblById, Order)
  }

	ensemblBySymbol1 <- NA
	ensemblBySymbol2 <- NA
  if (!.IsEmpty(geneSymbols)) {
		ensemblBySymbol1 <- .QueryEnsembl(inputIds = geneSymbols,
                             queryField = "external_gene_name",
                             dataset = dataset,
                             ensemblVersion = ensemblVersion,
														 ensemblMirror = ensemblMirror,
                             biomart = biomart)
		ensemblBySymbol1$GeneSymbol <- ensemblBySymbol1$external_gene_name
		ensemblBySymbol1 <- merge(combinedEnsembl, ensemblBySymbol1, by = 'GeneSymbol', all.x = T)
		ensemblBySymbol1 <- dplyr::arrange(ensemblBySymbol1, Order)

		ensemblBySymbol2 <- .QueryEnsembl(inputIds = geneSymbols,
														queryField = "hgnc_symbol",
														dataset = dataset,
														ensemblVersion = ensemblVersion,
														ensemblMirror = ensemblMirror,
														biomart = biomart)
		ensemblBySymbol2$GeneSymbol <- ensemblBySymbol2$hgnc_symbol
		ensemblBySymbol2 <- merge(combinedEnsembl, ensemblBySymbol2, by = 'GeneSymbol', all.x = T)
		ensemblBySymbol2 <- dplyr::arrange(ensemblBySymbol2, Order)
	}

	#Concat in preferential order, based on EnsemblId.  We can assume any resolved hit from Ensembl will have an Ensembl ID
	ret <- data.frame(EnsemblId = combinedEnsembl$EnsemblId, GeneSymbol = combinedEnsembl$GeneSymbol, ensembl_gene_id = NA, hgnc_symbol = NA, external_gene_name = NA, Order = 1:nrow(combinedEnsembl), stringsAsFactors=FALSE)
	ret <- .ConcatPreferentially(fieldToTest = 'ensembl_gene_id', datasets = list(ensemblById, ensemblBySymbol1, ensemblBySymbol2), baseDf = ret)
	ret <- ret[names(ret) != 'Order']

	return(ret)
}

.ConcatPreferentially <- function(fieldToTest, datasets, baseDf){
	datasets <- datasets[!is.na(datasets)]
	if (!('Order' %in% names(baseDf))) {
		baseDf$Order <- 1:nrow(baseDf)
	}

	ret <- baseDf[FALSE, TRUE]  #zero rows

	for (dataset in datasets) {
		toAppend <- dataset[!is.na(dataset[[fieldToTest]]) & !(dataset[[fieldToTest]] %in% ret[[fieldToTest]]),]
		ret <- rbind(ret, toAppend[names(ret)])
	}

	#now ensure all input terms represented, in order:
	ret <- rbind(ret, baseDf[!(baseDf$Order %in% ret$Order),])
	ret <- dplyr::arrange(ret, Order)

  return(ret)
}


# Translate a list of gene IDs to Ensembl IDs, based on a target field
.QueryEnsembl <- function(inputIds, queryField, biomart, dataset, ensemblVersion, ensemblMirror){
	print(paste0('Querying Ensembl using: ', queryField))
	ensembl <- biomaRt::useEnsembl(biomart = biomart,
		dataset = dataset,
		version = ensemblVersion,
		mirror = ensemblMirror
	)

	ensemblResults <- biomaRt::getBM(
		attributes = c('ensembl_gene_id', 'hgnc_symbol', 'external_gene_name'),
		filters = c(queryField),
		values = inputIds,
		mart = ensembl
	)

	#Drop duplicates.  Note: could consider group_concat on variables?
	ensemblResults <- ensemblResults %>% group_by_at(queryField) %>% mutate(total = dplyr::n())
	ensemblResults <- ensemblResults[ensemblResults$total == 1,]
	ensemblResults <- ensemblResults[names(ensemblResults) != 'total']

	print(paste0('Found ', sum(!is.na(ensemblResults$ensembl_gene_id)), ' of ', length(inputIds)))

	return(ensemblResults)
}


#' @title TranslateToStringDb
#' @param ensemblIds A vector of ensembl IDs
#' @param geneSymbols A vector of gene symbols
#' @param speciesId Species ID. see Stringdb reference for list of avialable species
#' @param replaceUnmatched Logical. If TRUE, removes NA's and replaces with inputs
#' @import STRINGdb
#' @importFrom dplyr %>% group_by summarize mutate_all coalesce
#' @importFrom BiocGenerics order
#' @importFrom stringr str_match
#' @importFrom stats setNames
TranslateToStringDb <- function(ensemblIds = NULL, geneSymbols = NULL, speciesId = 9606){
	.CheckGeneInputs(ensemblIds, geneSymbols)
  if (is.null(ensemblIds)) {
    ensemblIds <- NA
  }

  if (is.null(geneSymbols)) {
    geneSymbols <- NA
  }

  queryIds <- character()
	if (!.IsEmpty(ensemblIds)) {
    queryIds <- c(queryIds, ensemblIds)
  }

	if (!.IsEmpty(geneSymbols)) {
    queryIds <- c(queryIds, geneSymbols)
  }

  queryIds <- unique(queryIds)
  stringRes <- .QuerySTRINGdb(inputIds = queryIds, speciesId = speciesId)

	#Preferentially accept results by EnsemblId
	resultsBase <- data.frame(EnsemblId = ensemblIds, GeneSymbol = geneSymbols, stringsAsFactors = FALSE)
	resultsBase$Order <- 1:nrow(resultsBase)

	resultsById <- merge(resultsBase, stringRes, by.x = 'EnsemblId', all.x = T, by.y = 'InputTerm')
	resultsById <- resultsById[!is.na(resultsById$STRING_id),]

	resultsBySymbol <- merge(resultsBase, stringRes, by.x = 'GeneSymbol', all.x = T, by.y = 'InputTerm')
	resultsBySymbol <- resultsBySymbol[!(resultsBySymbol$Order %in% resultsById$Order),]

	results <- rbind(resultsById, resultsBySymbol)
	results <- rbind(results, resultsBase[!(resultsBase$Order %in% results$Order),])
	results <- dplyr::arrange(results, Order)
	results <- results[names(results) != 'Order']

	colnames(results)[colnames(results) == 'STRING_id'] <- 'STRING.id'

	return(results)
}



#' @title aliasTable
#' @param ensemblIds A vector of ensembl IDs
#' @param geneSymbols A vector of gene symbols
#' @param ensemblAttributes A vector of ensembl attributes
#' @param ensemblDataset Passed directly to biomaRt::useEnsembl
#' @param ensemblVersion Passed directly to biomaRt::useEnsembl
#' @param ensemblMirror Passed directly to biomaRt::useEnsembl
#' @param biomart Passed directly to biomaRt::useEnsembl
#' @param stringSpeciesId species ID. see Stringdb for list of available species
#' @param aliasPriorityOrder vector containig priority order of alias database. Must be UPPERCASE. Current databases: ENSEMBL, STRING
#' @importFrom biomaRt useEnsembl getBM
#' @import STRINGdb
#' @importFrom dplyr %>% group_by summarize mutate_all coalesce
#' @importFrom BiocGenerics order
#' @importFrom stringr str_match
#' @importFrom stats setNames
TranslateGeneNames <- function(ensemblIds = NULL, geneSymbols = NULL, ensemblDataset = "mmulatta_gene_ensembl", ensemblVersion = NULL, ensemblMirror = "uswest", biomart = 'ensembl',
                       stringSpeciesId = 9606,
											 useEnsembl = TRUE, useSTRINGdb = TRUE){

	.CheckGeneInputs(ensemblIds, geneSymbols)
	if (is.null(ensemblIds)) {
		ensemblIds <- NA
	}

	if (is.null(geneSymbols)) {
		geneSymbols <- NA
	}

	inputDf <- data.frame(EnsemblId = ensemblIds, GeneSymbol = geneSymbols, stringsAsFactors = FALSE)
	inputDf$Order <- 1:nrow(inputDf)

	if (useEnsembl) {
	  ret.ensembl <- TranslateToEnsembl(ensemblIds = ensemblIds,
                              geneSymbols = geneSymbols,
                              dataset = ensemblDataset,
                              ensemblVersion = ensemblVersion,
															ensemblMirror = ensemblMirror,
                              biomart = biomart)

		if (nrow(inputDf) != nrow(ret.ensembl)) {
			stop('Rows not equal for ensembl result')
		}

		ret.ensembl$Order <- 1:nrow(ret.ensembl)
		ret.ensembl <- ret.ensembl[!(names(ret.ensembl) %in% c('EnsemblId', 'GeneSymbol'))]
		inputDf <- merge(inputDf, ret.ensembl, by = 'Order', all.x = TRUE)
	}

	if (useSTRINGdb) {
  	ret.string <- TranslateToStringDb(ensemblIds = ensemblIds,
                              geneSymbols = geneSymbols,
                              speciesId = stringSpeciesId)

		if (nrow(inputDf) != nrow(ret.string)) {
			stop('Rows not equal for STRINGdb result')
		}

		ret.string$Order <- 1:nrow(ret.string)
		ret.string <- ret.string[!(names(ret.string) %in% c('EnsemblId', 'GeneSymbol'))]
		inputDf <- merge(inputDf, ret.string, by = 'Order', all.x = TRUE)
	}

	inputDf <- dplyr::arrange(inputDf, Order)
	inputDf <- inputDf[names(inputDf) != 'Order']

  return(inputDf)
}


