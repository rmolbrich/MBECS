# Sys.setenv('_R_CHECK_SYSTEM_CLOCK_' = 0)



# DEFINE CLASS ------------------------------------------------------------



#' Define MbecData-class
#'
#' An extension of phyloseq-class that contains the additional fields 'type', 'log' and
#' 'transformations' to accommodate MBECS functionality.
#' @keywords MBECS Class
#' @slot type User defined denominator for this data, e.g., 'RAW', 'normalized'.
#' @slot log log that will be filled by the other package functions
#' @slot input.obj either class phyloseq or a matrix of counts
#' @slot meta.obj dataframe of covariate variables
#' @slot required.col vector of strings that denote required variables in the covariate information
#' @slot tax_table taxonomic table from phyloseq as optional input
#' @slot phy_tree phylogenetic tree as optional input
#' @slot refseq reference sequences as optional input
#' @slot transformations list to be filled with the result of different correction methods
#' @import phyloseq
#' @export
#'
#' @examples
#' # use constructor with default parameters to create object from count-matrix
#' # and meta-data table.
#' mbec.obj <- MbecData(input.obj=datadummy$cnts, meta.obj = datadummy$meta)
MbecData <- setClass("MbecData", contains = "phyloseq",slots = list(type="character", log="character", transformations="list"))

#' Mbec-Data Constructor
#'
#' Constructor for the package class MbecData that takes a single input object of class phyloseq or a matrix of counts
#' and a data-frame of covariate variables for model-building.
#' @keywords MBECS Constructor
#' @param type string type that describes the data, e.g., raw, processed, ..
#' @param log log that will be filled by the other package functions
#' @param input.obj either class phyloseq or a matrix of counts
#' @param meta.obj A table with covariate information, whose row-names correspond to sample-IDs
#' @param required.col vector of strings that denote required variables in the covariate information
#' @param tax_table taxonomic table from phyloseq as optional input
#' @param phy_tree phylogenetic tree as optional input
#' @param refseq reference sequences as optional input
#' @param transformations list to be filled with the result of different correction methods
#' @return produces an R-object of type MbecData
#' @import phyloseq
#' @export
#'
#' @examples
#' # use constructor with default parameters to create object from count-matrix
#' # and meta-data table.
#' mbec.obj <- MbecData(input.obj=datadummy$cnts, meta.obj = datadummy$meta)
MbecData <- function(type=character(),
                     log=character(),
                     input.obj,
                     meta.obj=NULL,
                     required.col=c("group","batch"),
                     tax_table=NULL,
                     phy_tree=NULL,
                     refseq=NULL,
                     transformations=list()) {

  # if input is of class phyloseq - just include the type and done.
  if( is(input.obj, "phyloseq") ) {
    return( new("MbecData", type=type, log=log, input.obj) )
  } else {
    # if we got a matrix of counts - the meta object also needs to contain data

    ## Make all the necessary tests for matrix-like inputs and return counts
    if( !all(apply(input.obj, 2, is.numeric)) ) {
      stop("All the values in your input need to be numeric!", call. = FALSE)
    }
    ## check if meta is set and contains all the required columns
    if( is.null(meta.obj) | !all(eval(required.col) %in% (colnames(meta.obj))) ) {
      stop("You need to supply a meta-frame that contains the columns: ", paste(required.col, collapse=", "), call. = FALSE)
    }
    # check orientation of counts for phyloseq-constructor - meta.obj needs to be sxf anyway (if not then FU)
    if( dim(input.obj)[1] == length(meta.obj[,eval(required.col[1])]) ) {
      # taxa are columns
      return( new("MbecData", type=type, log=log, phyloseq::phyloseq(phyloseq::otu_table(input.obj, taxa_are_rows = FALSE),
                                                                     phyloseq::sample_data(meta.obj),
                                                                     phyloseq::tax_table(tax_table, errorIfNULL = FALSE),
                                                                     phyloseq::phy_tree(phy_tree, errorIfNULL = FALSE),
                                                                     phyloseq::refseq(refseq, errorIfNULL = FALSE))) )
    } else {
      return( new("MbecData", type=type, log=log, phyloseq::phyloseq(phyloseq::otu_table(input.obj, taxa_are_rows = TRUE),
                                                                     phyloseq::sample_data(meta.obj),
                                                                     phyloseq::tax_table(tax_table, errorIfNULL = FALSE),
                                                                     phyloseq::phy_tree(phy_tree, errorIfNULL = FALSE),
                                                                     phyloseq::refseq(refseq, errorIfNULL = FALSE))) )
    }
  }
}



# MbecData Setter ---------------------------------------------------------


#' Mbec-Data Setter
#'
#' This function either updates counts, type and log attributes (DEFAULT) or adds a matrix of
#' transformed counts to the input (update=FALSE).
#' @keywords MBECS Setter
#' @param input.obj MbecData object to change
#' @param new.cnts matrix-like object with same dimension as 'otu_table' in input.obj
#' @param log character to add to the log string
#' @param type character to replace type-attribute or as name tag for added count-matrix
#' @param update logical (TRUE) replace input fields (FALSE) add to the input.obj
#' @return Input object with updated attributes.
#' @export
#'
#' @examples
#' # This will replace the current abundance-matrix, append 'rbe' (for remove-batch-effect) to the
#' # log and change the type to 'rbe' to indicate current status.
#' MBEC.obj <- mbecSetData(input.obj=datadummy, new.cnts=datadummy$cnts,
#'     log='rbe',type='corrected', update=TRUE)
#'
#' # This will add the corrected data to the list of transformations named like the type argument.
#' # The 'log' attribute will be updated, but 'type' stays the same (because the main abundnce
#' # matrix wasn't replaced.
#' MBEC.obj <- mbecSetData(input.obj=datadummy, new.cnts=datadummy$cnts,
#'     log='rbe',type='corrected', update=FALSE)
mbecSetData <- function(input.obj, new.cnts=NULL, log=character(), type=character(), update=TRUE) {

  if( is.null(new.cnts) ) {
    message("Nothing to do - returning input unchanged.")
  }

  # UPDATE replaces values and appends log
  if( update ) {
    message("Update is: ", update,", replacing counts and setting type to '", type, "'.")
    # figure out orientation
    if( all(phyloseq::sample_names(input.obj) %in% colnames(new.cnts)) ) {
      # fxs orientation
      attr(input.obj, "otu_table") <- phyloseq::otu_table(new.cnts, taxa_are_rows = TRUE)

    } else if( all(phyloseq::sample_names(input.obj) %in% rownames(new.cnts)) ) {
      # sxf orientation
      attr(input.obj, "otu_table") <- phyloseq::otu_table(new.cnts, taxa_are_rows = FALSE)
    }

    # change type
    attr(input.obj, "type") <- type
    # append log
    tmp.log <- attr(input.obj, "log")
    message("old log is: ",tmp.log)
    tmp.log <- paste(tmp.log, log, sep = "|")
    attr(input.obj, "log") <- tmp.log
    message("new log is: ", tmp.log)

  } else {
    # no-UPDATE puts values in a list and does nothing to the log
    message("Update is: ", update,", adding counts to transformations-list.")

    # put new counts into the transformations list
    attr(input.obj, "transformations")[[eval(type)]] <- new.cnts

    # append log
    tmp.log <- attr(input.obj, "log")
    message("old log is: ",tmp.log)
    tmp.log <- paste(tmp.log, log, sep = "|")
    attr(input.obj, "log") <- tmp.log
    message("new log is: ", tmp.log)
  }
  return(input.obj)
}



# MbecData Getter ---------------------------------------------------------

#' Mbec-Data Getter
#'
#' This function extracts abundance matrix and meta-data in the chosen orientation from the input.
#'
#' The parameter 'orientation' determines if the output has features as columns (sxf) or if the
#' columns contain samples (fxs). This is mainly used to retrieve correctly oriented matrices for
#' the different analysis and correction functions.
#'
#' The parameter 'required.col' is a vector of column names (technically positions would work) from
#' metadata, that are required for the analysis at hand. The function actually only checks if
#' they are present in the data, but it will return the whole meta-frame.
#'
#' @keywords MBECS Getter
#' @param input.obj list(cnts, meta), phyloseq, MbecData object (correct orientation is handeled internally)
#' @param orientation, Select either 'fxs' or 'sxf' to retrieve features in rows or columns respectively
#' @param required.col Vector of column names that are required from the covariate-table.
#' @return A list that contains count-matrix (in chosen orientation) and meta-data table.
#' @export
#'
#' @examples
#' # This will return the abundance matrix with samples as rows and check for the presence of
#' # variables 'group' and 'batch' in the meta-data.
#' list.obj <- mbecGetData(input.obj=datadummy, orientation="sxf",
#'     required.col=c("group","batch"))
#'
#' # This will return the abundance matrix with samples as columns and check for the presence
#' # of variables 'group' and 'batch' in the meta-data.
#' list.obj <- mbecGetData(input.obj=datadummy, orientation="fxs",
#'     required.col=c("group","batch"))
setGeneric("mbecGetData", signature="input.obj",
           function(input.obj, orientation="fxs", required.col=NULL)
             standardGeneric("mbecGetData")
)

#IMPLEMENT function to use a vector as input
.mbecGetData <- function(input.obj, orientation="fxs", required.col=NULL) {

  # in general the input should be of class phyloseq or MbecData - so just get the fields to evaluate here
  tmp.meta <- data.frame(phyloseq::sample_data(input.obj, errorIfNULL = FALSE))

  # enable check of covariates if 'required.col' is supplied
  if( !is.null(required.col) ) {
    ## check if meta is set and contains all the required columns
    if( is.null(tmp.meta) | !all(eval(required.col) %in% (colnames(tmp.meta))) ) {
      stop("You need to supply a meta-frame that contains the columns: ", paste(required.col, collapse=", "), call. = FALSE)
    }
  }

  # 1. if everything is fine
  if( (!attr(phyloseq::otu_table(input.obj), "taxa_are_rows") & orientation == "sxf") || (attr(phyloseq::otu_table(input.obj), "taxa_are_rows") & orientation == "fxs") ) {
    tmp.cnts <- data.frame(phyloseq::otu_table(input.obj))

  } else if( (!attr(phyloseq::otu_table(input.obj), "taxa_are_rows") & orientation == "fxs") || (attr(phyloseq::otu_table(input.obj), "taxa_are_rows") & orientation == "sxf") ) {
    # if orientations do not fit - transpose count-matrix
    tmp.cnts <- data.frame(t(phyloseq::otu_table(input.obj)))

  } else {
    ### If this happens, sth. is severely broken.
    stop("Stop the shenanigans!")
  }

  # meta has always the same orientation
  tmp.meta <- data.frame(phyloseq::sample_data(input.obj))

  return(list(tmp.cnts, tmp.meta))

}


#' Mbec-Data Getter Phyloseq
#'
#' This function extracts abundance matrix and meta-data in the chosen orientation from input of
#' class phyloseq.
#'
#' The parameter 'orientation' determines if the output has features as columns (sxf) or if the
#' columns contain samples (fxs). This is mainly used to retrieve correctly oriented matrices for
#' the different analysis and correction functions.
#'
#' The parameter 'required.col' is a vector of column names (technically positions would work) from
#' metadata, that are required for the analysis at hand. The function actually only checks if
#' they are present in the data, but it will return the whole meta-frame.
#'
#' @keywords MBECS Getter
#' @param input.obj phyloseq-object
#' @param orientation, Select either 'fxs' or 'sxf' to retrieve features in rows or columns respectively
#' @param required.col Vector of column names that are required from the covariate-table.
#' @return A list that contains count-matrix (in chosen orientation) and meta-data table.
#' @export
#'
#' @examples
#' # This will return the abundance matrix with samples as rows and check for the presence of
#' # variables 'group' and 'batch' in the meta-data.
#' list.obj <- mbecGetData(input.obj=datadummy, orientation="sxf",
#'     required.col=c("group","batch"))
#'
#' # This will return the abundance matrix with samples as columns and check for the presence
#' # of variables 'group' and 'batch' in the meta-data.
#' list.obj <- mbecGetData(input.obj=datadummy, orientation="fxs",
#'     required.col=c("group","batch"))
setMethod("mbecGetData", "phyloseq",
          function(input.obj, orientation="fxs", required.col=NULL) {
            .mbecGetData(input.obj, orientation=orientation, required.col=required.col)
          }
)


#' Mbec-Data Getter MbecDdata
#'
#' This function extracts abundance matrix and meta-data in the chosen orientation from input of
#' class MbecData
#'
#' The parameter 'orientation' determines if the output has features as columns (sxf) or if the
#' columns contain samples (fxs). This is mainly used to retrieve correctly oriented matrices for
#' the different analysis and correction functions.
#'
#' The parameter 'required.col' is a vector of column names (technically positions would work) from
#' metadata, that are required for the analysis at hand. The function actually only checks if
#' they are present in the data, but it will return the whole meta-frame.
#'
#' @keywords MBECS Getter
#' @param input.obj MbecData-object
#' @param orientation, Select either 'fxs' or 'sxf' to retrieve features in rows or columns respectively
#' @param required.col Vector of column names that are required from the covariate-table.
#' @return A list that contains count-matrix (in chosen orientation) and meta-data table.
#' @export
#'
#' @examples
#' # This will return the abundance matrix with samples as rows and check for the presence of
#' # variables 'group' and 'batch' in the meta-data.
#' list.obj <- mbecGetData(input.obj=datadummy, orientation="sxf",
#'     required.col=c("group","batch"))
#'
#' # This will return the abundance matrix with samples as columns and check for the presence
#' # of variables 'group' and 'batch' in the meta-data.
#' list.obj <- mbecGetData(input.obj=datadummy, orientation="fxs",
#'     required.col=c("group","batch"))
setMethod("mbecGetData", "MbecData",
          function(input.obj, orientation="fxs", required.col=NULL) {
            .mbecGetData(input.obj, orientation=orientation, required.col=required.col)
          }
)


#' Mbec-Data Getter List
#'
#' This function extracts abundance matrix and meta-data in the chosen orientation from list input.
#' The implementation for lists basically performs some type checks on the input and then compares
#' col/row-names of abundance table and meta-data to infer orientation and return the user specified
#' direction
#'
#' The parameter 'orientation' determines if the output has features as columns (sxf) or if the
#' columns contain samples (fxs). This is mainly used to retrieve correctly oriented matrices for
#' the different analysis and correction functions.
#'
#' The parameter 'required.col' is a vector of column names (technically positions would work) from
#' meta-data, that are required for the analysis at hand. The function actually only checks if
#' they are present in the data, but it will return the whole meta-frame.
#'
#' @keywords MBECS Getter
#' @param input.obj phyloseq-object
#' @param orientation, Select either 'fxs' or 'sxf' to retrieve features in rows or columns respectively
#' @param required.col Vector of column names that are required from the covariate-table.
#' @return A list that contains count-matrix (in chosen orientation) and meta-data table.
#' @export
#'
#' @examples
#' # This will return the abundance matrix with samples as rows and check for the presence of
#' # variables 'group' and 'batch' in the meta-data.
#' list.obj <- mbecGetData(input.obj=datadummy, orientation="sxf",
#'     required.col=c("group","batch"))
#'
#' # This will return the abundance matrix with samples as columns and check for the presence
#' # of variables 'group' and 'batch' in the meta-data.
#' list.obj <- mbecGetData(input.obj=datadummy, orientation="fxs",
#'     required.col=c("group","batch"))
setMethod("mbecGetData", "list",
          function(input.obj, orientation="fxs", required.col=NULL) {
            message("It is a list!")
            ### for a list input with counts and meta-data, we need to do some more checks
            tmp.cnts <- input.obj[[1]]
            tmp.meta <- data.frame(input.obj[[2]])

            ## Make all the necessary tests for matrix-like inputs and return counts
            if( !all(apply(tmp.cnts, 2, is.numeric)) ) {
              stop("All the values in your input need to be numeric!", call. = FALSE)
            }

            # enable check of covariates if 'required.col' is supplied
            if( !is.null(required.col) ) {
              ## check if meta is set and contains all the required columns
              if( is.null(tmp.meta) | !all(eval(required.col) %in% (colnames(tmp.meta))) ) {
                stop("You need to supply a meta-frame that contains the columns: ", paste(required.col, collapse=", "), call. = FALSE)
              }
            }

            # figure out orientation - meta is supposed to have samplename in rows
            if( all(rownames(tmp.meta) %in% colnames(tmp.cnts)) ) {
              # fxs orientation - if it matches the 'orientation' flag just return, else transpose first
              if( orientation == "sxf" ) {
                tmp.cnts <- t(tmp.cnts)
              } # nothing else matters

            } else if( all(rownames(tmp.meta) %in% rownames(tmp.cnts)) ) {
              # sxf orientation - if it matches the 'orientation' flag just return, else transpose first
              if( orientation == "fxs" ) {
                tmp.cnts <- t(tmp.cnts)
              } # nothing else matters

            } else {
              ### If this happens, sth. is severely broken.
              stop("Stop the shenanigans! Either the dimensions of 'count' and 'meta' mtrices do not fit, or you lack col/row-names in your inputs.")
            }

            return(list(data.frame(tmp.cnts), tmp.meta))
          }
)



# mbecProcessInput --------------------------------------------------------


#' Mbec-Data Constructor Wrapper
#'
#' This function is a wrapper for the constructor of MbecData-objects. Given the parameter
#' 'required.col', the function will check if the required columns are present in the data and then
#' return it as an MbecData-object.
#'
#' The parameter 'required.col' is a vector of column names (technically positions would work) from
#' meta-data, that are required for the analysis at hand. The function actually only checks if
#' they are present in the data, but it will return the whole meta-frame.
#'
#' @keywords MBECS Constructor Wrapper
#' @param input.obj phyloseq-object
#' @param required.col Vector of column names that are required from the covariate-table.
#' @return An object of type MbecData that has been validated.
#' @export
#'
#' @examples
#' # This will check for the presence of variables 'group' and 'batch' in the meta-data and return
#' # an object of class 'MbecData'.
#' MbecData.obj <- mbecProcessInput(input.obj=datadummy,
#'     required.col=c("group","batch"))
setGeneric("mbecProcessInput", valueClass="MbecData", signature="input.obj",
           function(input.obj, required.col=NULL)
             standardGeneric("mbecProcessInput")
)

# the generic version is for MbecData type inputs
.mbecProcessInput <- function(input.obj, required.col=NULL) {

  # check sample_data if 'required.col' is not NULL
  # enable check of covariates if 'required.col' is supplied
  if( !is.null(required.col) ) {
    ## check if meta is set and contains all the required columns
    if( is.null(phyloseq::sample_data(input.obj, errorIfNULL = FALSE)) | !all(eval(required.col) %in% (phyloseq::sample_variables(input.obj))) ) {
      stop("You need to supply a meta-frame that contains the columns: ", paste(required.col, collapse=", "), call. = FALSE)
    }
  }
  # If validation succeeds, just return the input unchanged.
  return(input.obj)
}


#' Mbec-Data Constructor Wrapper for MbecData
#'
#' This function is a wrapper for the constructor of MbecData-objects. Given the parameter
#' 'required.col', the function will check if the required columns are present in the data and then
#' return it as an MbecData-object.
#'
#' The parameter 'required.col' is a vector of column names (technically positions would work) from
#' meta-data, that are required for the analysis at hand. The function actually only checks if
#' they are present in the data, but it will return the whole meta-frame.
#'
#' For existing MbecData objects this function serves as sanity check within correction and analysis
#' functions. Also kind of required to enable the input of different data classes at any point in
#' the MBECS pipeline.
#'
#' @keywords MBECS Constructor Wrapper MbecData
#' @param input.obj phyloseq-object
#' @param required.col Vector of column names that are required from the covariate-table.
#' @return An object of type MbecData that has been validated.
#' @export
#'
#' @examples
#' # This will check for the presence of variables 'group' and 'batch' in the meta-data and return
#' # an object of class 'MbecData'.
#' MbecData.obj <- mbecProcessInput(input.obj=datadummy,
#'     required.col=c("group","batch"))
setMethod("mbecProcessInput", "MbecData",
          function(input.obj, required.col=NULL) {
            .mbecProcessInput(input.obj, required.col=required.col)
          }
)


#' Mbec-Data Constructor Wrapper for Phyloseq
#'
#' This function is a wrapper for the constructor of MbecData-objects. Given the parameter
#' 'required.col', the function will check if the required columns are present in the data and then
#' return it as an MbecData-object.
#'
#' The parameter 'required.col' is a vector of column names (technically positions would work) from
#' meta-data, that are required for the analysis at hand. The function actually only checks if
#' they are present in the data, but it will return the whole meta-frame.
#'
#' @keywords MBECS Constructor Wrapper Phyloseq
#' @param input.obj phyloseq-object
#' @param required.col Vector of column names that are required from the covariate-table.
#' @return An object of type MbecData that has been validated.
#' @export
#'
#' @examples
#' # This will check for the presence of variables 'group' and 'batch' in the meta-data and return
#' # an object of class 'MbecData'.
#' MbecData.obj <- mbecProcessInput(input.obj=datadummy,
#'     required.col=c("group","batch"))
setMethod("mbecProcessInput", "phyloseq",
          function(input.obj, required.col=NULL) {

            # check sample_data if 'required.col' is not NULL
            # enable check of covariates if 'required.col' is supplied
            if( !is.null(required.col) ) {
              ## check if meta is set and contains all the required columns
              if( is.null(phyloseq::sample_data(input.obj, errorIfNULL = FALSE)) | !all(eval(required.col) %in% (phyloseq::sample_variables(input.obj))) ) {
                stop("You need to supply a meta-frame that contains the columns: ", paste(required.col, collapse=", "), call. = FALSE)
              }
            }

            # if everything checks out - create new MbecData object and return
            return.obj <- new("MbecData",
                              type="raw",
                              log="newobject",
                              transformations=list(),
                              input.obj)

            return(return.obj)
          }
)


#' Mbec-Data Constructor Wrapper for List Input
#'
#' This function is a wrapper for the constructor of MbecData-objects. Given the parameter
#' 'required.col', the function will check if the required columns are present in the data and then
#' return it as an MbecData-object.
#'
#' The parameter 'required.col' is a vector of column names (technically positions would work) from
#' meta-data, that are required for the analysis at hand. The function actually only checks if
#' they are present in the data, but it will return the whole meta-frame.
#'
#' @keywords MBECS Constructor Wrapper Phyloseq
#' @param input.obj a list that contains an abundance matrix and a data.frame of covaraite information
#' @param required.col Vector of column names that are required from the covariate-table.
#' @return An object of type MbecData that has been validated.
#' @export
#'
#' @examples
#' # This will check for the presence of variables 'group' and 'batch' in the meta-data and return
#' # an object of class 'MbecData'.
#' MbecData.obj <- mbecProcessInput(input.obj=datadummy,
#'     required.col=c("group","batch"))
setMethod("mbecProcessInput", "list",
          function(input.obj, required.col=NULL) {

            if( length(input.obj) != 2 ) {
              stop("Stop: Please provide the required input.")
            }

            # with numeric inputs, just call mbecGetData to handle all the testing and such
            input.obj <- mbecGetData(input.obj, orientation = "fxs", required.col = required.col)

            return.obj <- new("MbecData",
                              type="raw",
                              log="newobject",
                              transformations=list(),
                              phyloseq::phyloseq(phyloseq::otu_table(input.obj[[1]], taxa_are_rows = TRUE),
                                                 phyloseq::sample_data(input.obj[[2]])))

            return(return.obj)
          }
)

