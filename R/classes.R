##### CLASSES #####

#' CNV.anno class
#' @description Annotations required for CNV analysis are stored in this class.
#' @return \code{CNV.anno} class.
#' @details This class does not contain any sample data. Use \code{CNV.create_anno} to create.
#' @examples
#' # create object
#' anno <- CNV.create_anno()
#'
#' # general information
#' anno
#' show(anno)
#' @author Volker Hovestadt \email{conumee@@hovestadt.bio}
#' @export
setClass("CNV.anno", representation(date = "character", args = "list",
    genome = "data.frame", gap = "GRanges", probes = "GRanges", exclude = "GRanges",
    detail = "GRanges", bins = "GRanges"))

#' @rdname CNV.anno-class
#' @importFrom methods show
#' @param object \code{CNV.anno} object
#' @export
setMethod("show", "CNV.anno", function(object) {
    cat("CNV annotation object\n")
    cat("   created  : ", object@date, "\n", sep = "")
    cat("  @genome   : ", nrow(object@genome), " chromosomes\n", sep = "")
    cat("  @gap      : ", length(object@gap), " regions\n", sep = "")
    cat("  @probes   : ", length(object@probes), " probes\n", sep = "")
    cat("  @exclude  : ", length(object@exclude), " regions (overlapping ",
        length(findOverlaps(object@probes, object@exclude)), " probes)\n",
        sep = "")
    cat("  @detail   : ", length(object@detail), " regions (overlapping ",
        length(findOverlaps(object@probes, object@detail)), " probes)\n",
        sep = "")
    cat("  @bins     : ", length(object@bins), " bins (min/avg/max size: ",
        object@args$bin_minsize/1000, "/", suppressWarnings(round(mean(width(object@bins))/1000,
            1)), "/", object@args$bin_maxsize/1000, "kb, probes: ", object@args$bin_minprobes,
        "/", suppressWarnings(round(mean(values(object@bins)$probes), 1)),
        "/", max(values(object@bins)$probes), ")\n", sep = "")
})


#' CNV.data class
#' @description Intensities of one or multiple samples are stored in this class.
#' @return \code{CNV.data} class.
#' @details Use \code{CNV.load} to create.
#' @examples
#' # create object
#' library(minfiData)
#' data(MsetEx)
#'
#' d <- CNV.load(MsetEx)
#'
#' # general information
#' d
#' show(d)
#'
#' # show or replace sample names
#' names(d)
#' names(d) <- toupper(names(d))
#'
#' # subset samples
#' d[1:2]
#' @author Volker Hovestadt \email{conumee@@hovestadt.bio}
#' @export
setClass("CNV.data", representation(date = "character", intensity = "data.frame"))

#' @rdname CNV.data-class
#' @param object \code{CNV.data} object
setMethod("show", "CNV.data", function(object) {
    cat("CNV data object\n")
    cat("   created   : ", object@date, "\n", sep = "")
    if (length(object@intensity) == 0) {
        cat("  @intensity : unavailable, run CNV.load\n", sep = "")
    } else {
        cat("  @intensity : available (", ncol(object@intensity), " samples, ",
            nrow(object@intensity), " probes)\n", sep = "")
    }
})

#' @rdname CNV.data-class
#' @param x \code{CNV.data} object (defined by \code{Extract} generic).
#' @param i index. \code{logical}, \code{numeric} or \code{character}.
#' @export
setMethod("[", signature(x = "CNV.data"), function(x, i) {
    x@intensity <- x@intensity[, i, drop = FALSE]
    return(x)
})

#' @rdname CNV.data-class
setMethod("names", signature(x = "CNV.data"), function(x) {
    return(colnames(x@intensity))
})

#' @rdname CNV.data-class
#' @param value Replacement names.
setReplaceMethod("names", signature(x = "CNV.data"), function(x, value) {
    if (length(value) == ncol(x@intensity)) {
        colnames(x@intensity) <- value
    } else {
        stop("number of names does not fit number of samples.")
    }
    return(x)
})


#' CNV.analysis class
#' @description CNV analysis data of a single sample is stored in this class
#' @return \code{CNV.analysis} class.
#' @details Use \code{CNV.fit} to create. Modified by \code{CNV.bin}, \code{CNV.detail} and \code{CNV.segment}.
#' @examples
#' # prepare
#' library(minfiData)
#' data(MsetEx)
#' d <- CNV.load(MsetEx)
#' anno <- CNV.create_anno()
#'
#' # create object
#' x <- CNV.fit(query = d['GroupB_1'], ref = d[c('GroupA_1', 'GroupA_2', 'GroupA_3')], anno)
#'
#' # modify object
#' x <- CNV.bin(x)
#' x <- CNV.detail(x)
#' x <- CNV.segment(x)
#'
#' # general information
#' x
#' show(x)
#'
#' # coefficients of linear regression
#' coef(x)
#'
#' # show or replace sample name
#' names(x)
#' names(x) <- 'Sample 1'
#'
#' # output plots
#' CNV.genomeplot(x)
#' CNV.genomeplot(x, chr = 'chr6')
#' #CNV.detailplot(x, name = 'MYCN')
#' #CNV.detailplot_wrap(x)
#' CNV.write(x, what = 'segments')
#' @author Volker Hovestadt \email{conumee@@hovestadt.bio}
#' @export
setClass("CNV.analysis", representation(name = "character", date = "character",
    anno = "CNV.anno", fit = "list", bin = "list", detail = "list", seg = "list"))

#' @rdname CNV.analysis-class
#' @param object \code{CNV.analysis} object
setMethod("show", "CNV.analysis", function(object) {
    cat("CNV analysis object\n")
    cat("   created   : ", object@date, "\n", sep = "")
    cat("  @name      : ", object@name, "\n", sep = "")
    cat("  @anno      : ", nrow(object@anno@genome), " chromosomes, ",
        length(object@anno@probes), " probes, ", length(object@anno@bins),
        " bins\n", sep = "")
    if (length(object@fit) == 0) {
        cat("  @fit       : unavailable, run CNV.fit\n", sep = "")
    } else {
        cat("  @fit       : available (noise: ", round(object@fit$noise,
            3), ")\n", sep = "")
    }
    if (length(object@bin) == 0) {
        cat("  @bin       : unavailable, run CNV.bin\n", sep = "")
    } else {
        cat("  @bin       : available (shift: ", round(object@bin$shift,
            3), ")\n", sep = "")
    }
    if (length(object@detail) == 0) {
        cat("  @detail    : unavailable, run CNV.detail\n", sep = "")
    } else {
        cat("  @detail    : available (", length(object@detail$ratio),
            " regions)\n", sep = "")
    }
    if (length(object@seg) == 0) {
        cat("  @seg       : unavailable, run CNV.segment\n", sep = "")
    } else {
        cat("  @seg       : available (", nrow(object@seg$summary), " segments)\n",
            sep = "")
    }
})

#' @rdname CNV.analysis-class
#' @export
setMethod("names", signature(x = "CNV.analysis"), function(x) {
    x@name
})

#' @rdname CNV.analysis-class
#' @param x \code{CNV.analysis} object (defined by \code{show} generic).
#' @param value Replacement names.
#' @export
setReplaceMethod("names", signature(x = "CNV.analysis"), function(x, value) {
    if (length(value) == 1) {
        x@name <- value
    } else {
        stop("need exactly one sample name.")
    }
    return(x)
})

#' @rdname CNV.analysis-class
#' @importFrom stats coef
#' @export
setMethod("coef", signature(object = "CNV.analysis"), function(object) {
    object@fit$coef
})


#' CNV.summaryanalisys class
#' @description CNV summary analysis of grouped data of samples is stored in this class
#' @return \code{CNV.summaryanalysis} class.
#' @details Use \code{CNV.processsummary} to create.
#' @examples

#' # load RGset #
#' RGset<-readRDS('path/to/RGset.rds')
#' ### or use an Mset instead ###
#'
#' # process to Mset #
#' Mset<-preprocessIllumina(RGset)
#'
#' # load conumee package and load via CNV.load method #
#' library(conumee)
#' Data.data<- CNV.load(Mset)

#' # create annotation object; not complete and non working example below!! --> for reference and help read conumee vignette #
#' anno <- CNV.create_anno(array_type = array_type, detail_regions = detail_regions)
#'
#'
#' # fetch controls from your Data #
#' Data.controls<- grep( 'Control', names(Data.data))
#'
#' # read in phenosheet of samples; an identifier column mapping CNV data to phenosheet (named ID here ) should be specified #
#' test_pheno=read.csv(file ='../test_pheno.csv')
#' keep <- match(test_pheno$ID, names(Data.data))
#' Data.data_test<-Data.data[c(keep)]
#'
#' x<-CNV.processsummary(object=Data.data_test,
#'                                 pheno=test_pheno,
#'                                 labels='labels',
#'                                 interest_groups=NULL,
#'                                 identifier='ID',
#'                                 controls=c('Control'),
#'                                 anno= anno,
#'                                 summary_plots=TRUE,
#'                                 sample_plots=TRUE,
#'                                 chr = "all",
#'                                 chrX = TRUE,
#'                                 chrY = TRUE,
#'                                 centromere = TRUE,
#'                                 main = NULL,
#'                                 ylim = c(-1, 1),
#'                                 set_par = TRUE,
#'                                 save=TRUE,
#'                                 path="Path/to/savings")
#'
#' # general information
#' x
#' show(x)
#'
#' # show or replace sample name
#' names(x)
#' names(x) <- 'Sample 1'
#'
#' @author Samir Jabari \email{samir.jabari@@fau.de}
#' @export
setClass("CNV.summaryanalysis", representation(names = "character", date = "character", intensity_vals = "list",
cnv_seg_data = "list", cnvs = "list", gl_freqs= "list"
        ))

#' @rdname CNV.summaryanalysis-class
#' @param object \code{CNV.summaryanalysis} object
setMethod("show", "CNV.summaryanalysis", function(object) {
        cat("CNV summaryanalysis object\n")
        cat("   created   : ", object@date, "\n", sep = "")
        cat("  @names      : ", object@names, "\n", sep = "")

        if (length(object@intensity_vals) == 0 | is.null(object@intensity_vals)) {
            cat("  @intensity_vals       : unavailable, run CNV.processsummary\n", sep = "")
        } else {
            cat("  @intensity_vals       : available (groups: ", length(object@intensity_vals), ")\n", sep = "")
        }

        if (length(object@cnv_seg_data) == 0 | is.null(object@cnv_seg_data)) {
            cat("  @cnv_seg_data       : unavailable, run CNV.processsummary\n", sep = "")
        } else {
            cat("  @cnv_seg_data       : available (groups: ", length(object@cnv_seg_data), ")\n", sep = "")
        }

        if (length(object@cnvs) == 0 | is.null(object@cnvs)) {
            cat("  @cnvs       : unavailable, run CNV.processsummary\n", sep = "")
        } else {
            cat("  @cnvs       : available (groups: ", length(object@cnvs), ")\n", sep = "")
        }

        if (length(object@gl_freqs) == 0 | is.null(object@gl_freqs)) {
            cat("  @gl_freqs       : unavailable, run CNV.processsummary\n", sep = "")
        } else {
            cat("  @gl_freqs       : available (groups: ", length(object@gl_freqs), ")\n", sep = "")
        }

    })

#' @rdname CNV.summaryanalysis-class
#' @export
setMethod("names", signature(x = "CNV.summaryanalysis"), function(x) {
        x@names
    })

#' @rdname CNV.summaryanalysis-class
#' @param x \code{CNV.summaryanalysis} object (defined by \code{show} generic).
#' @param value Replacement names.
#' @export
setReplaceMethod("names", signature(x = "CNV.summaryanalysis"), function(x, value) {
        if (length(value) == 1) {
            x@name <- value
        } else {
            stop("need exactly one sample name.")
        }
        return(x)
    })
