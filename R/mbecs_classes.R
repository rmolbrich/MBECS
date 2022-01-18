# DEFINE CLASS ------------------------------------------------------------


#' Define MbecData-class
#'
#' An extension of phyloseq-class that contains the additional attributes 'tss',
#' 'clr', 'corrections' and 'assessments' to enable the MBECS functionality.
#'
#' @keywords MBECS Class
#' @slot otu_table Class phyloseq::otu_table, (usually sparse) matrix of
#' abundance values.
#' @slot sample_data Dataframe of covariate variables.
#' @slot tax_table Taxonomic table from phyloseq as optional input
#' @slot phy_tree Phylogenetic tree as optional input
#' @slot refseq Reference sequences as optional input
#' @slot assessments A list for the results of batch effect assessment
#' algorithms (BEAA) that produce p-values for all features.
#' @slot corrections A list for the results of batch effect correction
#' algorithms (BECA) that produce adjusted abundance matrices.
#' @slot tss Total-sum-squared feature abundance matrix.
#' @slot clr Cumulative log-ratio transformed feature abundance matrix.
#' @import phyloseq
#' @export
#'
#' @examples
#' # use constructor with default parameters to create object from count-matrix
#' # and meta-data table.
#' mbec.obj <- MbecData(cnt_table=dummy.list$cnts, meta_data = dummy.list$meta)
MbecData <- setClass("MbecData", contains = "phyloseq",
                     slots = list(assessments="list", corrections="list",
                                  tss="matrix",       clr="matrix"))


#' Mbec-Data Constructor
#'
#' Constructor for the package class MbecData. Minimum input is an abundance
#' matrix for the argument 'cnt_table' and any type of frame that contains
#' columns of covariate information. The argument 'cnt_table' requires col/row-
#' names that correspond to features and samples. The correct orientation will
#' be handled internally. The argument 'meta_data' requires row-names that
#' correspond to samples. Although it is an exported function, the user should
#' utilize the function 'mbecProcessInput()' for safe initialization of an
#' MbecData-object from phyloseq or list(counts, metadata) inputs.
#'
#' Additional (OPTIONAL) arguments are 'tax_table', 'phy_tree' and 'ref_seq'
#' from phyloseq-objects.
#'
#' The (OPTIONAL) arguments 'tss' and 'clr' are feature
#' abundance matrices that should contain total-sum-scaled or cumulative
#' log-ratio transformed counts respectively. They should however be calculated
#' by the package-function 'mbecTransform()'.
#'
#' The lists for Assessments and Corrections will be initialized empty and
#' should only be accessed via the available Get/Set-functions.
#'
#' @keywords MBECS Constructor
#' @param cnt_table either class phyloseq or a matrix of counts
#' @param meta_data A table with covariate information, whose row-names correspond to sample-IDs
#' @param tax_table taxonomic table from phyloseq as optional input
#' @param phy_tree phylogenetic tree as optional input
#' @param refseq reference sequences as optional input
#' @param assessments A list for the results of BEAAs.
#' @param corrections A list for the results of BECAs.
#' @param tss Total-sum-squared features matrix.
#' @param clr Cumulative log-ratio transformed feature matrix.
#' @return produces an R-object of type MbecData
#' @import phyloseq
#' @export
#'
#' @examples
#' # use constructor with default parameters to create object from count-matrix
#' # and meta-data table.
#' mbec.obj <- MbecData(cnt_table=dummy.list$cnts, meta_data = dummy.list$meta)
MbecData <- function( cnt_table=NULL,
                      meta_data=NULL,
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
    if( all(meta_data$sID %in% colnames(tss)) ) {
      # fxs orientation do nothing

    } else if( all(meta_data$sID %in% rownames(tss)) ) {
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
    if( all(meta_data$sID %in% colnames(clr)) ) {
      # fxs orientation do nothing

    } else if( all(meta_data$sID %in% rownames(clr)) ) {
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
      return( new("MbecData",
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
      return( new("MbecData",
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
#' Sets and/or replaces selected feature abundance matrix and handles correct
#' orientation. The argument type determines which slot to access, i.e. the base
#' matrices for un-transformed counts "otu", total sum-scaled counts "tss",
#' cumulative log-ratio transformed counts "clr" and batch effect corrected
#' counts "cor" and assessment vectors "ass". The later two additionally require
#' the use of the argument 'label' that specifies the name within the respective
#' lists of corrections and assessments.
#'
#' @keywords MBECS Setter
#' @param input.obj MbecData object to work on.
#' @param new.cnts A matrix-like object with same dimension as 'otu_table' in
#' input.obj.
#' @param type Specify which type of data to add, by using one of 'ass'
#' (Assessement), 'cor' (Correction), 'clr' (Cumulative Log-Ratio) or 'tss'
#' (Total Scaled-Sum).
#' @param label For types 'ass' and 'cor' this sets the name within the lists.
#' @return Input object with updated attributes.
#' @export
#'
#' @examples
#' # This will fill the 'tss' slot with the supplied matrix.
#' MBEC.obj <- mbecSetData(input.obj=dummy.mbec, new.cnts=dummy.list$cnts,
#'     type='tss')
#'
#' # This will put the given matrix into the list of corrected counts under the
#' # name "nameOfMethod".
#' MBEC.obj <- mbecSetData(input.obj=dummy.mbec, new.cnts=dummy.list$cnts,
#'     type='cor', label="nameOfMethod")
setGeneric("mbecSetData", signature="input.obj",
           function(input.obj, new.cnts=NULL, type=c("otu","ass","cor","clr","tss"),
                    label=character())
             standardGeneric("mbecSetData")
)


#' Mbec-Data Setter
#'
#' Sets and/or replaces selected feature abundance matrix and handles correct
#' orientation. The argument type determines which slot to access, i.e. the base
#' matrices for un-transformed counts "otu", total sum-scaled counts "tss",
#' cumulative log-ratio transformed counts "clr" and batch effect corrected
#' counts "cor" and assessment vectors "ass". The later two additionally require
#' the use of the argument 'label' that specifies the name within the respective
#' lists of corrections and assessments.
#'
#' @keywords MBECS Setter
#' @param input.obj MbecData object to work on.
#' @param new.cnts A matrix-like object with same dimension as 'otu_table' in
#' input.obj.
#' @param type Specify which type of data to add, by using one of 'ass'
#' (Assessement), 'cor' (Correction), 'clr' (Cumulative Log-Ratio) or 'tss'
#' (Total Scaled-Sum).
#' @param label For types 'ass' and 'cor' this sets the name within the lists.
#' @return Input object with updated attributes.
#' @export
#'
#' @examples
#' # This will fill the 'tss' slot with the supplied matrix.
#' MBEC.obj <- mbecSetData(input.obj=dummy.mbec, new.cnts=dummy.list$cnts,
#'     type='tss')
#'
#' # This will put the given matrix into the list of corrected counts under the
#' # name "nameOfMethod".
#' MBEC.obj <- mbecSetData(input.obj=dummy.mbec, new.cnts=dummy.list$cnts,
#'     type='cor', label="nameOfMethod")
setMethod("mbecSetData", "MbecData",
          function(input.obj, new.cnts=NULL,
                   type=c("otu","ass","cor","clr","tss"), label=character()) {
          .mbecSetData(input.obj, new.cnts=new.cnts,
                       type=type, label=label)
          }
)


#' Mbec-Data Setter
#'
#' Sets and/or replaces selected feature abundance matrix and handles correct
#' orientation. The argument type determines which slot to access, i.e. the base
#' matrices for un-transformed counts "otu", total sum-scaled counts "tss",
#' cumulative log-ratio transformed counts "clr" and batch effect corrected
#' counts "cor" and assessment vectors "ass". The later two additionally require
#' the use of the argument 'label' that specifies the name within the respective
#' lists of corrections and assessments.
#'
#' @keywords MBECS Setter
#' @param input.obj MbecData object to work on.
#' @param new.cnts A matrix-like object with same dimension as 'otu_table' in
#' input.obj.
#' @param type Specify which type of data to add, by using one of 'ass'
#' (Assessement), 'cor' (Correction), 'clr' (Cumulative Log-Ratio) or 'tss'
#' (Total Scaled-Sum).
#' @param label For types 'ass' and 'cor' this sets the name within the lists.
#' @return Input object with updated attributes.
#' @export
#'
#' @examples
#' # This will fill the 'tss' slot with the supplied matrix.
#' MBEC.obj <- mbecSetData(input.obj=dummy.mbec, new.cnts=dummy.list$cnts,
#'     type='tss')
#'
#' # This will put the given matrix into the list of corrected counts under the
#' # name "nameOfMethod".
#' MBEC.obj <- mbecSetData(input.obj=dummy.mbec, new.cnts=dummy.list$cnts,
#'     type='cor', label="nameOfMethod")
.mbecSetData <- function(input.obj, new.cnts=NULL,
                        type=c("otu","ass","cor","clr","tss"), label=character()) {

  if( is.null(new.cnts) ) {
    warning("Nothing to do here. Returning unchanged input.")
    return(input.obj)
  }

  type <- match.arg(type, choices = c("otu","ass","cor","clr","tss"))


  if( type != "ass" ) {
    if( !all(apply(new.cnts, 1, is.numeric)) ) {
      stop("All the values in your count-table need to be numeric!", call. = FALSE)
    } else if( is.null(colnames(new.cnts)) && is.null(rownames(new.cnts)) ) {
      stop("No col/row-names found! Not able to match to sample names.")
    }
  } else {
    if( !all(apply(new.cnts, 1, is.numeric)) ) {
      stop("All the values in your assessment vector need to be numeric!", call. = FALSE)
    } else if( is.null(names(new.cnts)) ) {
      stop("No col/row-names found! Not able to match to features.")
    }

  }

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
      label <- paste("item",(length(input.obj@assessments)+1), sep="_")
    input.obj@assessments[[eval(label)]] <- new.cnts
  } else if( type == "cor" ) {
    if( length(label) == 0 )
      label <- paste("item",(length(input.obj@corrections)+1), sep="_")
    input.obj@corrections[[eval(label)]] <- new.cnts
  }

  return(input.obj)
}



# MbecData Getter ---------------------------------------------------------

#' Mbec-Data Getter
#'
#' This function extracts abundance matrix and meta-data in the chosen
#' orientation from the input.
#'
#' The parameter 'orientation' determines if the output has features as columns
#' (sxf) or if the columns contain samples (fxs). This is mainly used to
#' retrieve correctly oriented matrices for the different analysis and
#' correction functions.
#'
#' The parameter 'required.col' is a vector of column names (technically
#' positions would work) in the  metadata, that are required for the analysis
#' at hand. The function actually only checks if they are present in the data,
#' but it will return the whole meta-frame.
#'
#' The argument type determines which slot to access, i.e. the base
#' matrices for un-transformed counts "otu", total sum-scaled counts "tss",
#' cumulative log-ratio transformed counts "clr" and batch effect corrected
#' counts "cor" and assessment vectors "ass". The later two additionally require
#' the use of the argument 'label' that specifies the name within the respective
#' lists of corrections and assessments.
#'
#' @keywords MBECS Getter
#' @param input.obj MbecData object
#' @param orientation, Select either 'fxs' or 'sxf' to retrieve features in
#' rows or columns respectively.
#' @param required.col Vector of column names that are required from the
#' covariate-table.
#' @param type Specify which type of data to add, by using one of 'ass'
#' (Assessement), 'cor' (Correction), 'clr' (Cumulative Log-Ratio) or 'tss'
#' (Total Scaled-Sum).
#' @param label For types 'ass' and 'cor' this specifies the name within the
#' lists.
#' @return A list that contains count-matrix (in chosen orientation) and
#' meta-data table.
#' @export
#'
#' @examples
#' # This will return the un-transformed (OTU) abundance matrix with features as
#' # columns and it will test if the columns "group" and "batch" are present in
#' # the meta-data table.
#' list.obj <- mbecGetData(input.obj=dummy.mbec, orientation="sxf",
#'     required.col=c("group","batch"), type="otu")
#'
#' # This will return the clr-transformed abundance matrix with features as
#' # rows and it will test if the columns "group" and "batch" are present in
#' # the meta-data table.
#' list.obj <- mbecGetData(input.obj=dummy.mbec, orientation="fxs",
#'     required.col=c("group","batch"), type="clr")
setGeneric("mbecGetData", signature="input.obj",
           function(input.obj, orientation="fxs", required.col=NULL,
                    type=c("otu","ass","cor","clr","tss"), label=character())
             standardGeneric("mbecGetData")
)


#' Mbec-Data Getter
#'
#' This function extracts abundance matrix and meta-data in the chosen
#' orientation from the input.
#'
#' The parameter 'orientation' determines if the output has features as columns
#' (sxf) or if the columns contain samples (fxs). This is mainly used to
#' retrieve correctly oriented matrices for the different analysis and
#' correction functions.
#'
#' The parameter 'required.col' is a vector of column names (technically
#' positions would work) in the  metadata, that are required for the analysis
#' at hand. The function actually only checks if they are present in the data,
#' but it will return the whole meta-frame.
#'
#' The argument type determines which slot to access, i.e. the base
#' matrices for un-transformed counts "otu", total sum-scaled counts "tss",
#' cumulative log-ratio transformed counts "clr" and batch effect corrected
#' counts "cor" and assessment vectors "ass". The later two additionally require
#' the use of the argument 'label' that specifies the name within the respective
#' lists of corrections and assessments.
#'
#' @keywords MBECS Getter
#' @param input.obj MbecData object
#' @param orientation, Select either 'fxs' or 'sxf' to retrieve features in
#' rows or columns respectively.
#' @param required.col Vector of column names that are required from the
#' covariate-table.
#' @param type Specify which type of data to add, by using one of 'ass'
#' (Assessement), 'cor' (Correction), 'clr' (Cumulative Log-Ratio) or 'tss'
#' (Total Scaled-Sum).
#' @param label For types 'ass' and 'cor' this specifies the name within the
#' lists.
#' @return A list that contains count-matrix (in chosen orientation) and
#' meta-data table.
#' @export
#'
#' @examples
#' # This will return the un-transformed (OTU) abundance matrix with features as
#' # columns and it will test if the columns "group" and "batch" are present in
#' # the meta-data table.
#' list.obj <- mbecGetData(input.obj=dummy.mbec, orientation="sxf",
#'     required.col=c("group","batch"), type="otu")
#'
#' # This will return the clr-transformed abundance matrix with features as
#' # rows and it will test if the columns "group" and "batch" are present in
#' # the meta-data table.
#' list.obj <- mbecGetData(input.obj=dummy.mbec, orientation="fxs",
#'     required.col=c("group","batch"), type="clr")
.mbecGetData <- function(input.obj, orientation="fxs", required.col=NULL,
                         type=c("otu","ass","cor","clr","tss"), label=character()) {

  type <- match.arg(type, choices = c("otu","ass","cor","clr","tss"))
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
    if( length(label) == 0 || !(label %in% names(input.obj@assessments)) ) {
      message("The given matrix (label parameter) is not part of this object.
              Returning the complete assessments list.")
      tmp.cnts <- input.obj@assessments
    } else {
      tmp.cnts <- data.frame(input.obj@assessments[[eval(label)]], check.names = FALSE)
    }
  } else if( type == "cor" ) {
    if( length(label) == 0 || !(label %in% names(input.obj@corrections)) ) {
      message("The given matrix (label parameter) is not part of this object.
              Returning the complete corrections list.")
      tmp.cnts <- input.obj@corrections
    } else {
      tmp.cnts <- data.frame(input.obj@corrections[[eval(label)]], check.names = FALSE)
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
    # This will just return the otu object.
  }

  return(list(tmp.cnts, tmp.meta))
}


#' Mbec-Data Getter
#'
#' This function extracts abundance matrix and meta-data in the chosen
#' orientation from the input.
#'
#' The parameter 'orientation' determines if the output has features as columns
#' (sxf) or if the columns contain samples (fxs). This is mainly used to
#' retrieve correctly oriented matrices for the different analysis and
#' correction functions.
#'
#' The parameter 'required.col' is a vector of column names (technically
#' positions would work) in the  metadata, that are required for the analysis
#' at hand. The function actually only checks if they are present in the data,
#' but it will return the whole meta-frame.
#'
#' The argument type determines which slot to access, i.e. the base
#' matrices for un-transformed counts "otu", total sum-scaled counts "tss",
#' cumulative log-ratio transformed counts "clr" and batch effect corrected
#' counts "cor" and assessment vectors "ass". The later two additionally require
#' the use of the argument 'label' that specifies the name within the respective
#' lists of corrections and assessments.
#'
#' @keywords MBECS Getter
#' @param input.obj MbecData object
#' @param orientation, Select either 'fxs' or 'sxf' to retrieve features in
#' rows or columns respectively.
#' @param required.col Vector of column names that are required from the
#' covariate-table.
#' @param type Specify which type of data to add, by using one of 'ass'
#' (Assessement), 'cor' (Correction), 'clr' (Cumulative Log-Ratio) or 'tss'
#' (Total Scaled-Sum).
#' @param label For types 'ass' and 'cor' this specifies the name within the
#' lists.
#' @return A list that contains count-matrix (in chosen orientation) and
#' meta-data table.
#' @export
#'
#' @examples
#' # This will return the un-transformed (OTU) abundance matrix with features as
#' # columns and it will test if the columns "group" and "batch" are present in
#' # the meta-data table.
#' list.obj <- mbecGetData(input.obj=dummy.mbec, orientation="sxf",
#'     required.col=c("group","batch"), type="otu")
#'
#' # This will return the clr-transformed abundance matrix with features as
#' # rows and it will test if the columns "group" and "batch" are present in
#' # the meta-data table.
#' list.obj <- mbecGetData(input.obj=dummy.mbec, orientation="fxs",
#'     required.col=c("group","batch"), type="clr")
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
#' This function is a wrapper for the constructor of MbecData-objects from
#' phyloseq objects and lists of counts and sample data.
#'
#' The (OPTIONAL) argument 'required.col' is a vector of column-names that will
#' enable a sanity test for the presence in the meta-data table. Which is also
#' the second use-case for objects that are already of class MbecData.
#'
#' @keywords MBECS Constructor Wrapper
#' @param input.obj One of MbecData, phyloseq or list(counts, meta-data).
#' @param required.col Vector of column names that need to be present in the
#' meta-data table.
#' @return An object of type MbecData that has been validated.
#' @export
#'
#' @examples
#' # This will check for the presence of variables 'group' and 'batch' in the
#' # meta-data and return an object of class 'MbecData'.
#' MbecData.obj <- mbecProcessInput(input.obj=dummy.mbec,
#'     required.col=c("group","batch"))
setGeneric("mbecProcessInput", valueClass="MbecData", signature="input.obj",
           function(input.obj, required.col=NULL)
             standardGeneric("mbecProcessInput")
)


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


#' Mbec-Data Constructor Wrapper
#'
#' This function is a wrapper for the constructor of MbecData-objects from
#' phyloseq objects and lists of counts and sample data.
#'
#' The (OPTIONAL) argument 'required.col' is a vector of column-names that will
#' enable a sanity test for the presence in the meta-data table. Which is also
#' the second use-case for objects that are already of class MbecData.
#'
#' @keywords MBECS Constructor Wrapper
#' @param input.obj One of MbecData, phyloseq or list(counts, meta-data).
#' @param required.col Vector of column names that need to be present in the
#' meta-data table.
#' @return An object of type MbecData that has been validated.
#' @export
#'
#' @examples
#' # This will check for the presence of variables 'group' and 'batch' in the
#' # meta-data and return an object of class 'MbecData'.
#' MbecData.obj <- mbecProcessInput(input.obj=dummy.mbec,
#'     required.col=c("group","batch"))
setMethod("mbecProcessInput", "MbecData",
          function(input.obj, required.col=NULL) {
            .mbecProcessInput(input.obj, required.col=required.col)
          }
)


#' Mbec-Data Constructor Wrapper
#'
#' This function is a wrapper for the constructor of MbecData-objects from
#' phyloseq objects and lists of counts and sample data.
#'
#' The (OPTIONAL) argument 'required.col' is a vector of column-names that will
#' enable a sanity test for the presence in the meta-data table. Which is also
#' the second use-case for objects that are already of class MbecData.
#'
#' @keywords MBECS Constructor Wrapper
#' @param input.obj One of MbecData, phyloseq or list(counts, meta-data).
#' @param required.col Vector of column names that need to be present in the
#' meta-data table.
#' @return An object of type MbecData that has been validated.
#' @export
#'
#' @examples
#' # This will check for the presence of variables 'group' and 'batch' in the
#' # meta-data and return an object of class 'MbecData'.
#' MbecData.obj <- mbecProcessInput(input.obj=dummy.mbec,
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
                                   tax_table = phyloseq::tax_table(input.obj, errorIfNULL = FALSE),
                                   phy_tree = phyloseq::phy_tree(input.obj, errorIfNULL = FALSE),
                                   refseq = phyloseq::refseq(input.obj, errorIfNULL = FALSE))

            return(return.obj)
          }
)


#' Mbec-Data Constructor Wrapper
#'
#' This function is a wrapper for the constructor of MbecData-objects from
#' phyloseq objects and lists of counts and sample data.
#'
#' The (OPTIONAL) argument 'required.col' is a vector of column-names that will
#' enable a sanity test for the presence in the meta-data table. Which is also
#' the second use-case for objects that are already of class MbecData.
#'
#' @keywords MBECS Constructor Wrapper
#' @param input.obj One of MbecData, phyloseq or list(counts, meta-data).
#' @param required.col Vector of column names that need to be present in the
#' meta-data table.
#' @return An object of type MbecData that has been validated.
#' @export
#'
#' @examples
#' # This will check for the presence of variables 'group' and 'batch' in the
#' # meta-data and return an object of class 'MbecData'.
#' MbecData.obj <- mbecProcessInput(input.obj=dummy.mbec,
#'     required.col=c("group","batch"))
setMethod("mbecProcessInput", "list",
          function(input.obj, required.col=NULL) {

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
                                   meta_data = input.obj[[2]])

            return(return.obj)
          }
)



# mbecGetPhyloseq ---------------------------------------------------------


#' Return Phyloseq after correction
#'
#' This function extracts the abundance table of choice and returns a phyloseq
#' object for downstream analyses.
#'
#' The argument type determines which slot to access, i.e. the base
#' matrices for un-transformed counts "otu", total sum-scaled counts "tss",
#' cumulative log-ratio transformed counts "clr" and batch effect corrected
#' counts "cor". The later additionally requires the use of the argument 'label'
#' that specifies the name within the list of corrected matrices.
#'
#' @keywords MBECS Phyloseq Getter
#' @param input.obj MbecData object
#' @param type Specify which type of data to add, by using one of 'cor'
#' (Correction), 'clr' (Cumulative Log-Ratio) or 'tss' (Total Scaled-Sum).
#' @param label For type 'cor' this specifies the name within the
#' list.
#' @return A phyloseq object that contains the chosen abundance table as
#' otu_table.
#' @export
#'
#' @examples
#' # This will return a phyloseq object that contains the clr-transformed
#' # abundances as otu_table
#' ps.clr.obj <- mbecGetPhyloseq(input.obj=dummy.mbec, type="clr")
setGeneric("mbecGetPhyloseq", signature="input.obj",
           function(input.obj, type=c("otu","cor","clr","tss"),
                    label=character())
             standardGeneric("mbecGetPhyloseq")
)



#' Return Phyloseq after correction
#'
#' This function extracts the abundance table of choice and returns a phyloseq
#' object for downstream analyses.
#'
#' The argument type determines which slot to access, i.e. the base
#' matrices for un-transformed counts "otu", total sum-scaled counts "tss",
#' cumulative log-ratio transformed counts "clr" and batch effect corrected
#' counts "cor". The later additionally requires the use of the argument 'label'
#' that specifies the name within the list of corrected matrices.
#'
#' @keywords MBECS Phyloseq Getter
#' @param input.obj MbecData object
#' @param type Specify which type of data to add, by using one of 'cor'
#' (Correction), 'clr' (Cumulative Log-Ratio) or 'tss' (Total Scaled-Sum).
#' @param label For type 'cor' this specifies the name within the
#' list.
#' @return A phyloseq object that contains the chosen abundance table as
#' otu_table.
#' @export
#'
#' @examples
#' # This will return a phyloseq object that contains the clr-transformed
#' # abundances as otu_table
#' ps.clr.obj <- mbecGetPhyloseq(input.obj=dummy.mbec, type="clr")
.mbecGetPhyloseq <- function(input.obj, type=c("otu","cor","clr","tss"), label=character()) {

  # match correct argument
  type <- match.arg(type, choices = c("otu","cor","clr","tss"))

  if( type == "cor" && length(label) == 0 ) {
    warning("For type 'cor' you need to specifiy the label-parameter.")
    return(NULL)
  }
  # extract table in fxs orientation
  table2use <- mbecGetData(input.obj, orientation = "fxs", type=eval(type), label=eval(label))[[1]]
  # create return object of class phyloseq
  ret.ps <- phyloseq::phyloseq(phyloseq::otu_table(table2use , taxa_are_rows = TRUE),
                               input.obj@sam_data,
                               input.obj@tax_table,
                               input.obj@phy_tree,
                               input.obj@refseq)
  # and return
  return(ret.ps)
}



#' Return Phyloseq after correction
#'
#' This function extracts the abundance table of choice and returns a phyloseq
#' object for downstream analyses.
#'
#' The argument type determines which slot to access, i.e. the base
#' matrices for un-transformed counts "otu", total sum-scaled counts "tss",
#' cumulative log-ratio transformed counts "clr" and batch effect corrected
#' counts "cor". The later additionally requires the use of the argument 'label'
#' that specifies the name within the list of corrected matrices.
#'
#' @keywords MBECS Phyloseq Getter
#' @param input.obj MbecData object
#' @param type Specify which type of data to add, by using one of 'cor'
#' (Correction), 'clr' (Cumulative Log-Ratio) or 'tss' (Total Scaled-Sum).
#' @param label For type 'cor' this specifies the name within the
#' list.
#' @return A phyloseq object that contains the chosen abundance table as
#' otu_table.
#' @export
#'
#' @examples
#' # This will return a phyloseq object that contains the clr-transformed
#' # abundances as otu_table
#' ps.clr.obj <- mbecGetPhyloseq(input.obj=dummy.mbec, type="clr")
setMethod("mbecGetPhyloseq", "MbecData",
          function(input.obj, type=c("otu","cor","clr","tss"),
                   label=character()) {
            .mbecGetPhyloseq(input.obj, type=type, label=label)
          }
)

