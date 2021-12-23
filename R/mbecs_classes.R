# Sys.setenv('_R_CHECK_SYSTEM_CLOCK_' = 0)

# DEFINE CLASS ------------------------------------------------------------


#' Define MbecData-class
#'
#' An extension of phyloseq-class that contains the additional fields 'type', 'log' and
#' 'transformations' to accommodate MBECS functionality.
#' @keywords MBECS Class
#' @slot type User defined denominator for this data, e.g., 'RAW', 'normalized'.
#' @slot log log that will be filled by the other package functions
#' @slot otu_table Class phyloseq::otu_table, (usually sparse) matrix of
#' abundance values.
#' @slot sample_data Dataframe of covariate variables.
#' @slot tax_table taxonomic table from phyloseq as optional input
#' @slot phy_tree phylogenetic tree as optional input
#' @slot refseq reference sequences as optional input
#' @slot assessments A list for the results of BEAAs.
#' @slot corrections A list for the results of BECAs.
#' @slot tss Total-sum-squared features matrix.
#' @slot clr Cumulative log-ratio transformed feature matrix.
#' @import phyloseq
#' @export
#'
#' @examples
#' # use constructor with default parameters to create object from count-matrix
#' # and meta-data table.
#' mbec.obj <- MbecData(cnt_table=datadummy$cnts, meta_data = datadummy$meta)
MbecData <- setClass("MbecData", contains = "phyloseq",
                     slots = list(type="character",   log="character",
                                  assessments="list", corrections="list",
                                  tss="matrix",       clr="matrix"))


#' Mbec-Data Constructor
#'
#' Constructor for the package class MbecData that takes a single input object of class phyloseq or a matrix of counts
#' and a data-frame of covariate variables for model-building.
#' @keywords MBECS Constructor
#' @param type string type that describes the data, e.g., raw, processed, ..
#' @param log log that will be filled by the other package functions
#' @param cnt_table either class phyloseq or a matrix of counts
#' @param meta_data A table with covariate information, whose row-names correspond to sample-IDs
#' @param tax_table taxonomic table from phyloseq as optional input
#' @param phy_tree phylogenetic tree as optional input
#' @param refseq reference sequences as optional input
#' @slot assessments A list for the results of BEAAs.
#' @slot corrections A list for the results of BECAs.
#' @slot tss Total-sum-squared features matrix.
#' @slot clr Cumulative log-ratio transformed feature matrix.
#' @return produces an R-object of type MbecData
#' @import phyloseq
#' @export
#'
#' @examples
#' # use constructor with default parameters to create object from count-matrix
#' # and meta-data table.
#' mbec.obj <- MbecData(cnt_table=datadummy$cnts, meta_data = datadummy$meta)
MbecData <- function( cnt_table=NULL,
                      meta_data=NULL,
                      type=character(),
                      log=character(),
                      tax_table=NULL,
                      phy_tree=NULL,
                      refseq=NULL,
                      assessments=list(),
                      corrections=list(),
                      tss=NULL,
                      clr=NULL) {

    ## Make all the necessary tests for matrix-like inputs and return counts
    if( !all(apply(cnt_table, 2, is.numeric)) ) {
      stop("All the values in your count-table need to be numeric!", call. = FALSE)
    }
    ## check if meta is set and contains all the required columns
    if( is.null(meta_data) ) {
      stop("You need to supply a meta-frame.", call. = FALSE)
    }
    ## check for 'sID' column and set it as row-names - or create 'sID' column from row-names
    if( any(colnames(meta_data) %in% "sID") ) {
      # 'sID' exists --> test if it is equal to rownames - else set as rownames
      if( !all(meta_data[,"sID"] %in% rownames(meta_data)) ) {
        message("sID column does not match with rownames, changing rownames to match sID column now.")
      } # else is ok


    } else {
      # check for rownames and basically make the taxa_are_rows_check here already
      message("No 'sID' column present, creating from rownames now.")
      meta_data$sID <- rownames(meta_data)
    }

  # IF initialized with 'tss' or 'clr' this will ensure 'fxs' orientation
  # DEBUG: tss and clr init needs some work but this will do for now.. I think

  if( !is.null(tss) ) {
    ## 1. figure out orientation and adjust accordingly
    if( all(phyloseq::sample_names(input.obj) %in% colnames(tss)) ) {
      # fxs orientation do nothing

    } else if( all(phyloseq::sample_names(input.obj) %in% rownames(tss)) ) {
      # change from sxf to fxs orientation
      tss <- t(tss)

      } else {
      warning("sth. is fishy with tss rotation")
    }
  } else {
    tss <- matrix(0)
  }

  if( !is.null(clr) ) {
    ## 1. figure out orientation and adjust accordingly
    if( all(phyloseq::sample_names(input.obj) %in% colnames(clr)) ) {
      # fxs orientation do nothing

    } else if( all(phyloseq::sample_names(input.obj) %in% rownames(clr)) ) {
      # change from sxf to fxs orientation
      clr <- t(clr)

    } else {
      warning("sth. is fishy with clr rotation")
    }
  } else {
    clr <- matrix(0)
  }


    # check orientation of counts for phyloseq-constructor - meta_data needs to be sxf anyway (if not then FU)
    if( dim(cnt_table)[1] == length(meta_data[,"sID"]) ) {
      # taxa are columns
      return( new("MbecData", type=type, log=log,
                  phyloseq::phyloseq(phyloseq::otu_table(cnt_table, taxa_are_rows = FALSE),
                                                                     phyloseq::sample_data(meta_data),
                                                                     phyloseq::tax_table(tax_table, errorIfNULL = FALSE),
                                                                     phyloseq::phy_tree(phy_tree, errorIfNULL = FALSE),
                                                                     phyloseq::refseq(refseq, errorIfNULL = FALSE)),
              assessments=assessments,
              corrections=corrections,
              tss=tss,
              clr=clr))
    } else {
      return( new("MbecData", type=type, log=log,
                  phyloseq::phyloseq(phyloseq::otu_table(cnt_table, taxa_are_rows = TRUE),
                                                                     phyloseq::sample_data(meta_data),
                                                                     phyloseq::tax_table(tax_table, errorIfNULL = FALSE),
                                                                     phyloseq::phy_tree(phy_tree, errorIfNULL = FALSE),
                                                                     phyloseq::refseq(refseq, errorIfNULL = FALSE)),
              assessments=assessments,
              corrections=corrections,
              tss=tss,
              clr=clr))
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
#' @param type Specify which type of data to add, by using one of 'ass' (Assessement), 'cor' (Correction), 'clr' (Cumulative Log-Ratio Transformation) or 'tss' (Total Scaled-Sum).
#' @param label For types 'ass' and 'cor' this sets the name within the list.
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
setGeneric("mbecSetData", signature="input.obj",
           function(input.obj, new.cnts=NULL, type=c("otu","ass","cor","clr","tss"),
                    label=character())
             standardGeneric("mbecSetData")
)


setMethod("mbecSetData", "MbecData",
          function(input.obj, new.cnts=NULL,
                   type=c("otu","ass","cor","clr","tss"), label=character()) {
          .mbecSetData(input.obj, new.cnts=new.cnts,
                       type=type, label=label)
          }
)


.mbecSetData <- function(input.obj, new.cnts=NULL,
                        type=c("otu","ass","cor","clr","tss"), label=character()) {

  if( is.null(new.cnts) ) {
    stop("Nothing to do here.")
  } else if( !all(apply(new.cnts, 2, is.numeric)) ) {
    stop("All the values in your count-table need to be numeric!", call. = FALSE)
  } else if( is.null(colnames(new.cnts)) && is.null(rownames(bla)) ) {
    stop("No col/row-names found! Not able to match to sample names.")
  }

  type <- match.arg(type, choices = c("otu","ass","cor","clr","tss"))

  # 1. get sample names
  s.names <- phyloseq::sample_names(input.obj)
  # 2. figure out orientation and change to 'fxs' if necessary
  if( all(s.names %in% colnames(new.cnts)) ) {
    # fxs orientation
    # do nothing

  } else if( all(s.names %in% rownames(new.cnts)) ) {
    # sxf orientation - rotate matrix
    new.cnts <- t(new.cnts)
  }

  if( type == "otu" ) {
    input.obj@otu_table <- phyloseq::otu_table(new.cnts, taxa_are_rows = TRUE)
  } else if( type == "clr" ) {
    input.obj@clr <- new.cnts
  } else if( type == "tss" ) {
    message("Set tss-transformed counts.")
    input.obj@tss <- new.cnts
  } else if( type == "ass" ) {
    if( length(label) == 0 )
      label <- paste("item",(length(input.obj@ass)+1), sep="_")
    input.obj@ass[[eval(label)]] <- new.cnts
  } else if( type == "cor" ) {
    if( length(label) == 0 )
      label <- paste("item",(length(input.obj@cor)+1), sep="_")
    input.obj@cor[[eval(label)]] <- new.cnts
  }

  return(input.obj)
}



# MbecData Getter ---------------------------------------------------------

#' Mbec-Data Getter
#'
#' This mow function extracts abundance matrix and meta-data in the chosen orientation from the input.
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
#' @param transformation Default is 'NULL', otherwise give name or index of
#' matrix in transformation-list attribute for extraction.
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
           function(input.obj, orientation="fxs", required.col=NULL,
                    type=c("otu","ass","cor","clr","tss"), label=character())
             standardGeneric("mbecGetData")
)

#IMPLEMENT function to use a vector as input
.mbecGetData <- function(input.obj, orientation="fxs", required.col=NULL,
                         type=c("otu","ass","cor","clr","tss"), label=character()) {

  type <- match.arg(type, choices = c("otu","ass","cor","clr","tss"))
  message("Selected type is: ", type)
  # in general the input should be of class phyloseq or MbecData - so just get the fields to evaluate here
  tmp.meta <- data.frame(phyloseq::sample_data(input.obj, errorIfNULL = FALSE), check.names = FALSE)

  # enable check of covariates if 'required.col' is supplied
  if( !is.null(required.col) ) {
    ## check if meta is set and contains all the required columns
    if( is.null(tmp.meta) || !all(eval(required.col) %in% (colnames(tmp.meta))) ) {
      stop("You need to supply a meta-frame that contains the columns: ", paste(required.col, collapse=", "), call. = FALSE)
    }
  }

  if( type == "otu" ) {
    tar <- attr(phyloseq::otu_table(input.obj), "taxa_are_rows")
    # 1. if everything is fine
    if( (!tar & orientation == "sxf") || (tar & orientation == "fxs") ) {
      tmp.cnts <- data.frame(as(phyloseq::otu_table(input.obj),"matrix"), check.names = FALSE)
    } else if( (!tar & orientation == "fxs") || (tar & orientation == "sxf") ) {
      # if orientations do not fit - transpose count-matrix
      tmp.cnts <- data.frame(as(t(phyloseq::otu_table(input.obj)),"matrix"), check.names = FALSE)
    }
  } else if( type == "clr" ) {
    tmp.cnts <- data.frame(input.obj@clr, check.names = FALSE)
  } else if( type == "tss" ) {
    tmp.cnts <- data.frame(input.obj@tss, check.names = FALSE)
  } else if( type == "ass" ) {
    if( length(label) == 0 || !(label %in% names(input.obj@ass)) ) {
      message("The given matrix (label parameter) is not part of this object.
              Returning the complete assessments list.")
      tmp.cnts <- input.obj@ass
    } else {
      tmp.cnts <- data.frame(input.obj@ass[[eval(label)]], check.names = FALSE)
    }
  } else if( type == "cor" ) {
    if( length(label) == 0 || !(label %in% names(input.obj@cor)) ) {
      message("The given matrix (label parameter) is not part of this object.
              Returning the complete corrections list.")
      tmp.cnts <- input.obj@cor
    } else {
      tmp.cnts <- data.frame(input.obj@cor[[eval(label)]], check.names = FALSE)
    }
  }

  ## 1. figure out orientation and adjust accordingly
  if( !is(tmp.cnts, "list") && (type != "otu") ) {
    if( (all(phyloseq::sample_names(input.obj) %in% colnames(tmp.cnts))
         && orientation == "sxf") ||
        (all(phyloseq::sample_names(input.obj) %in% rownames(tmp.cnts))
         && orientation == "fxs") ) {
      tmp.cnts <- data.frame(t(tmp.cnts), check.names = FALSE)
    }
  } else {
    # ToDo: this is not correct - FixMe: I work but the message is flawed
    message("Returning a list-object. Correct orientation will not be checked.")
  }

  return(list(tmp.cnts, tmp.meta))
}


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
#' @param transformation Default is 'NULL', otherwise give name or index of
#' matrix in transformation-list attribute for extraction.
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
          function(input.obj, orientation="fxs", required.col=NULL,
                   type=c("otu","ass","cor","clr","tss"), label=character()) {
            .mbecGetData(input.obj, orientation=orientation, required.col=required.col,
                         type=type, label=label)
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

            message("DEBUG: input for mbecProcessInput is of type phyloseq.")
            # check sample_data if 'required.col' is not NULL
            # enable check of covariates if 'required.col' is supplied
            if( !is.null(required.col) ) {
              ## check if meta is set and contains all the required columns
              if( is.null(phyloseq::sample_data(input.obj, errorIfNULL = FALSE)) | !all(eval(required.col) %in% (phyloseq::sample_variables(input.obj))) ) {
                stop("You need to supply a meta-frame that contains the columns: ", paste(required.col, collapse=", "), call. = FALSE)
              }
            }

            return.obj <- MbecData(cnt_table = phyloseq::otu_table(input.obj, taxa_are_rows = TRUE),
                                   meta_data = phyloseq::sample_data(input.obj),
                                   type="raw",
                                   log="newobject",
                                   tax_table = phyloseq::tax_table(input.obj, errorIfNULL = FALSE),
                                   phy_tree = phyloseq::phy_tree(input.obj, errorIfNULL = FALSE),
                                   refseq = phyloseq::refseq(input.obj, errorIfNULL = FALSE))

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

            message("DEBUG: input for mbecProcessInput is of type list")

            if( length(input.obj) != 2 ) {
              stop("Stop: Please provide an abundance-table as first element
                   and meta-data as second element of the list.")
            }

            # enable check of covariates if 'required.col' is supplied
            if( !is.null(required.col) ) {
              ## check if meta is set and contains all the required columns
              if( is.null(input.obj[[2]]) | !all(eval(required.col) %in% (colnames(input.obj[[2]]))) ) {
                stop("You need to supply a meta-frame that contains the columns: ", paste(required.col, collapse=", "), call. = FALSE)
              }
            }

            return.obj <- MbecData(cnt_table = input.obj[[1]],
                                   meta_data = input.obj[[2]],
                                   type="raw",
                                   log="newobject")

            return(return.obj)
          }
)

