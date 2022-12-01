##### PROCESSING methods #####

#' CNV.fit
#' @description Normalize query sample intensities by fitting intensities to reference set using a linear regression model.
#' @param query \code{CNV.data} object of query sample (single sample).
#' @param ref \code{CNV.data} object of reference set.
#' @param anno \code{CNV.anno} object. Use \code{CNV.create_anno} do create.
#' @param name character. Optional parameter to set query sample name.
#' @param intercept logical. Should intercept be considered? Defaults to \code{TRUE}.
#' @param ... Additional parameters (\code{CNV.fit} generic, currently not used).
#' @return \code{CNV.analysis} object.
#' @details The log2 ratio of query intensities versus a linear combination of reference set intensities that best reflects query intensities is calculated (as determined by linear regression). The annotations provided to \code{CNV.fit} are saved within the returned \code{CNV.analysis} object and used for subsequent analysis steps.
#' @examples
#' # prepare
#' library(minfiData)
#' data(MsetEx)
#' d <- CNV.load(MsetEx)
#' data(detail_regions)
#' anno <- CNV.create_anno(detail_regions = detail_regions)
#'
#' # create object
#' x <- CNV.fit(query = d['GroupB_1'], ref = d[c('GroupA_1', 'GroupA_2', 'GroupA_3')], anno)
#'
#' # modify object
#' #x <- CNV.bin(x)
#' #x <- CNV.detail(x)
#' #x <- CNV.segment(x)
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
#' @author Volker Hovestadt \email{conumee@@hovestadt.bio}
#' @export
setGeneric("CNV.fit", function(query, ref, anno, ...) {
    standardGeneric("CNV.fit")
})

#' @rdname CNV.fit
setMethod("CNV.fit", signature(query = "CNV.data", ref = "CNV.data", anno = "CNV.anno"),
    function(query, ref, anno, name = NULL, intercept = TRUE) {
        if (ncol(query@intensity) == 0)
            stop("query intensities unavailable, run CNV.load")
        if (ncol(ref@intensity) == 0)
            stop("reference set intensities unavailable, run CNV.load")

        if (ncol(query@intensity) != 1)
            stop("query contains more than one sample.")
        if (ncol(ref@intensity) == 1)
            warning("reference set contains only a single sample. use more samples for better results.")

        p <- names(anno@probes)  # ordered by location
        if (!all(is.element(p, rownames(query@intensity))))
            stop("query intensities not given for all probes.")
        if (!all(is.element(p, rownames(ref@intensity))))
            stop("reference set intensities not given for all probes.")

        object <- new("CNV.analysis")
        object@date <- date()
        object@fit$args <- list(intercept = intercept)

        if (!is.null(name)) {
            names(object) <- name
        } else {
            names(object) <- colnames(query@intensity)
        }
        object@anno <- anno

        r <- cor(query@intensity[p, ], ref@intensity[p, ])[1, ] < 0.99
        if (any(!r)) message("query sample seems to also be in the reference set. not used for fit.")
        if (intercept) {
            ref.fit <- lm(y ~ ., data = data.frame(y = query@intensity[p,
                1], X = ref@intensity[p, r]))
        } else {
            ref.fit <- lm(y ~ . - 1, data = data.frame(y = query@intensity[p,
                1], X = ref@intensity[p, r]))
        }
        object@fit$coef <- ref.fit$coefficients

        ref.predict <- predict(ref.fit)
        ref.predict[ref.predict < 1] <- 1

        object@fit$ratio <- log2(query@intensity[p, 1]/ref.predict[p])
        object@fit$noise <- sqrt(mean((object@fit$ratio[-1] - object@fit$ratio[-length(object@fit$ratio)])^2,
            na.rm = TRUE))

        return(object)
    })


#' CNV.bin
#' @description Combine single probe intensitiy values into predefined bins.
#' @param object \code{CNV.analysis} object.
#' @param ... Additional parameters (\code{CNV.bin} generic, currently not used).
#' @return \code{CNV.analysis} object.
#' @details The median intensity per bin is calculated. Bins are defined using \code{CNV.create_anno}. A value by which all probe and bin intensity values are shifted in subsequent analysis steps is calculated by minimizing the median absolute deviation from all bins to zero (ideally shifting the copy-neutral state to 0).
#' @examples
#' # prepare
#' library(minfiData)
#' data(MsetEx)
#' d <- CNV.load(MsetEx)
#' data(detail_regions)
#' anno <- CNV.create_anno(detail_regions = detail_regions)
#'
#' # create object
#' x <- CNV.fit(query = d['GroupB_1'], ref = d[c('GroupA_1', 'GroupA_2', 'GroupA_3')], anno)
#'
#' # modify object
#' x <- CNV.bin(x)
#' #x <- CNV.detail(x)
#' #x <- CNV.segment(x)
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
#' @author Volker Hovestadt \email{conumee@@hovestadt.bio}
#' @export
setGeneric("CNV.bin", function(object, ...) {
    standardGeneric("CNV.bin")
})

#' @rdname CNV.bin
setMethod("CNV.bin", signature(object = "CNV.analysis"), function(object) {
    if (length(object@fit) == 0)
        stop("fit unavailable, run CNV.fit")

    o1 <- as.matrix(findOverlaps(query = object@anno@bins, subject = object@anno@probes))
    o2 <- data.frame(bin = names(object@anno@bins)[o1[, "queryHits"]],
        probe = names(object@anno@probes)[o1[, "subjectHits"]], stringsAsFactors = FALSE)

    object@bin$ratio <- sapply(split(object@fit$ratio[o2[, "probe"]], o2[,
        "bin"]), median, na.rm = TRUE)[names(object@anno@bins)]
    object@bin$shift <- optim(0, function(s) median(abs(object@bin$ratio -
        s), na.rm = TRUE), method = "Brent", lower = -100, upper = 100)$par

    return(object)
})


#' CNV.detail
#' @description Combine single probe values within detail regions.
#' @param object \code{CNV.analysis} object.
#' @param ... Additional parameters (\code{CNV.detail} generic, currently not used).
#' @return \code{CNV.analysis} object.
#' @details The median intensity per detail region is calculated. Detail regions are defined using \code{CNV.create_anno(detail_bed=)}
#' @examples
#' # prepare
#' library(minfiData)
#' data(MsetEx)
#' d <- CNV.load(MsetEx)
#' data(detail_regions)
#' anno <- CNV.create_anno(detail_regions = detail_regions)
#'
#' # create object
#' x <- CNV.fit(query = d['GroupB_1'], ref = d[c('GroupA_1', 'GroupA_2', 'GroupA_3')], anno)
#'
#' # modify object
#' x <- CNV.bin(x)
#' x <- CNV.detail(x)
#' #x <- CNV.segment(x)
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
#' @author Volker Hovestadt \email{conumee@@hovestadt.bio}
#' @export
setGeneric("CNV.detail", function(object, ...) {
    standardGeneric("CNV.detail")
})

#' @rdname CNV.detail
setMethod("CNV.detail", signature(object = "CNV.analysis"), function(object) {
    if (length(object@fit) == 0)
        stop("fit unavailable, run CNV.fit")
    # if(length(object@bin) == 0) stop('bin unavailable, run CNV.bin')

    if (length(object@anno@detail) == 0) {
        message("no detail regions provided, define using CNV.create_anno")
    } else {
        d1 <- as.matrix(findOverlaps(query = object@anno@detail, subject = object@anno@probes))
        d2 <- data.frame(detail = values(object@anno@detail)$name[d1[,
            "queryHits"]], probe = names(object@anno@probes[d1[, "subjectHits"]]),
            stringsAsFactors = FALSE)

        object@detail$ratio <- sapply(split(object@fit$ratio[d2[, "probe"]],
            d2[, "detail"]), median, na.rm = TRUE)[values(object@anno@detail)$name]
        object@detail$probes <- table(d2[, 1])[values(object@anno@detail)$name]
    }
    return(object)
})


#' @import DNAcopy
NULL

#' CNV.segment
#' @description Segment bin values (wrapper of \code{DNAcopy} package).
#' @param object \code{CNV.analysis} object.
#' @param alpha See details. Defaults to 0.001.
#' @param nperm See details. Defaults to 50000.
#' @param min.width See details. Defaults to 5.
#' @param undo.splits See details. Defaults to 'sdundo'.
#' @param undo.SD See details. Defaults to 2.2.
#' @param verbose See details. Defaults to 0.
#' @param ... Additional parameters supplied to the \code{segment} method of the \code{DNAcopy} package.
#' @return \code{CNV.analysis} object.
#' @details This method is a wrapper of the CNA, segment, segments.summary and segments.p methods of the DNAcopy package. Please refer to the respective man pages for more detailed information. The default parameters of \code{CNV.segment} override some of the default parameters of segment and are optimized for 450k data CNV analysis.
#' @examples
#' # prepare
#' library(minfiData)
#' data(MsetEx)
#' d <- CNV.load(MsetEx)
#' data(detail_regions)
#' anno <- CNV.create_anno(detail_regions = detail_regions)
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
#' @author Volker Hovestadt \email{conumee@@hovestadt.bio}
#' @export
setGeneric("CNV.segment", function(object, ...) {
    standardGeneric("CNV.segment")
})

#' @rdname CNV.segment
setMethod("CNV.segment", signature(object = "CNV.analysis"), function(object,
    alpha = 0.001, nperm = 50000, min.width = 5, undo.splits = "sdundo",
    undo.SD = 2.2, verbose = 0, ...) {
    # if(length(object@fit) == 0) stop('fit unavailable, run CNV.fit')
    if (length(object@bin) == 0)
        stop("bin unavailable, run CNV.bin")
    # if(length(object@detail) == 0) stop('bin unavailable, run
    # CNV.detail')

    a1 <- formals()
    a2 <- as.list(match.call())[-1]
    object@seg$args <- as.list(sapply(setdiff(unique(names(c(a1, a2))),
        c("object", "verbose")), function(an) if (is.element(an, names(a2)))
        a2[[an]] else a1[[an]], simplify = FALSE))

    x1 <- DNAcopy::CNA(genomdat = object@bin$ratio[names(object@anno@bins)],
        chrom = as.vector(seqnames(object@anno@bins)), maploc = values(object@anno@bins)$midpoint,
        data.type = "logratio", sampleid = "sampleid")
    x2 <- DNAcopy::segment(x = x1, verbose = verbose, min.width = min.width,
        nperm = nperm, alpha = alpha, undo.splits = undo.splits, undo.SD = undo.SD,
        ...)
    object@seg$summary <- DNAcopy::segments.summary(x2)
    object@seg$summary$chrom <- as.vector(object@seg$summary$chrom)  # DNAcopy will factor chrom names. is there another way?
    object@seg$p <- DNAcopy::segments.p(x2)
    object@seg$p$chrom <- as.vector(object@seg$p$chrom)

    return(object)
})


#' CNV.processsummary
#' @description Create a CNV.summaryanalysis object containing grouped cnv analysis data with option of summaryplotting or single sample analysis and plots.
#' @param object \code{CNV.data} object.
#' @param pheno=NULL, data.frame phenosheet containing info about samples with columns as seen in other parameters
#' @param labels=NULL, charcter. labels column in pheno sheet
#' @param interest_groups=NULL, character specify the labels/groups you are interested in analysing (needs to be column name of phenosheet)
#' @param identifier=NULL, character. (phenosheet-)column mapping CNV data (from names()) to this phenosheet column
#' @param controls=NULL, character. labelname of controls in labels column
#' @param anno=NULL,  CNV.anno object. annotation object from CNV.create_anno
#' @param cutoff=0.4, numeric. threshold value to call gain or loss of glFrequency of DNAcopy package
#' @param summary_plots=TRUE, logical. whether to perform summary plots
#' @param sample_plots=TRUE, logical. whether to perform single sample plots (by now all single sample methods will be performed from the conumee package)
#' @param summaryanalysis CNV.summaryanalysis object as returned from CNV.processsummary function i.e.
#' @param intensity_vals. intensityvalues for the groups of interest specified as matrices as lists.
#' @param cnv_seg_data CNV.segment sample data as list grouped by groups of interest.
#' @param gl_freqs gl_frequency objects of the DNAcopy package for selected groups as lists.
#' @param chr = "all", character vector. Which chromomsomes to plot. Defaults to 'all'.
#' @param chrX = TRUE, logical. Plot values for chrX? Defaults to TRUE. Set CNV.create_anno(chrXY =FALSE)  if chrX and Y should not be included at all.
#' @param chrY = TRUE,  logical. Plot values for chrY? Defaults to TRUE
#' @param centromere = TRUE, logical. Show dashed lines at centromeres? Defaults to TRUE.
#' @param main = NULL, character. Title of the plot. Defaults to interest groups names.
#' @param ylim = c(-1, 1) numeric vector. The y limits of the plot. Defaults to c(-1, 1).
#' @param set_par logical. Use recommended graphical parameters for \code{oma} and \code{mar}? Defaults to \code{TRUE}. Original parameters are restored afterwards.
#' @param save= TRUE, logical. whether to save to path or not
#' @param path= NULL character. path to savings without trailing / !
#' @return \code{CNV.summaryanalysis} object.
#' @details This method generates a CNV.summaryanalyis object needed to generate summary plots with CNV.summaryplot of specified interest groups. On the y-axis gains and losses are plotted according to group frequency. See parameters for more information.
#' @examples
#' # see CNV.summaryanalisys for an example
#' @author Samir Jabari\email{samir.jabari@@fau.de}
#' @export

setGeneric("CNV.processsummary", function(object, ...) {
    standardGeneric("CNV.processsummary")
})
setMethod("CNV.processsummary", signature(object="CNV.data"), function(object=NULL, # CNV.data object
                                pheno=NULL, # phenosheet
                                labels=NULL,
                                interest_groups=NULL,
                                identifier=NULL,
                                controls=NULL,
                                anno=NULL,
                                cutoff=0.4,
                                summary_plots=TRUE,
                                sample_plots=TRUE,
                                chr = "all",
                                chrX = TRUE,
                                chrY = TRUE,
                                centromere = TRUE,
                                main = NULL,
                                ylim = c(-1, 1),
                                set_par = TRUE,
                                save=FALSE,
                                path=NULL){


if (summary_plots || sample_plots ){
    if (save){
      if (is.null(path)){
              stop('You need to specify a path' )
              }
      else {
      dir.create(path, recursive=TRUE, showWarnings = FALSE)
      }
    }
}


if (!class(object)=='CNV.data'){
     stop("Input must be of type 'CNV.data' ")
                               }

if (is.null(controls)){
    stop("Controls must be specified")
                       }

if (is.null(identifier)){
    stop("An identifier which connects the CNV.data and phenosheet data must be specified as column identifier from the phenosheet")
      }


if (is.null(labels)){
    stop("A 'labels' column of your phenosheet must be specified")
      }

if (is.null(anno)){
          stop("Please create and provide an annotation object using conumee's 'CNV.create_anno' method")
                    }


cnv_list<-list()
colnames_list<-list()
intensityvals_list<-list()

control<-pheno[grep(paste(controls,collapse='|'), pheno[[labels]]),]

if (is.null(interest_groups)){
    message('You did not specify a group of interest; using all groups available from the labels option specified!')

    interest_groups <- c(unique(pheno[grep(paste(controls,collapse='|'), pheno[[labels]], invert = TRUE),][labels]))[[1]]


     }


groups<-pheno[grep(paste(interest_groups,collapse='|'), pheno[[labels]]),]


keep <- match(control[[identifier]], names(object))
object.controls<-object[c(keep),]

keep <- match(groups[[identifier]], names(object) )

object.groups<-object[c(keep),]


for (i in seq_along(names(object.groups))){

        print (i)
        fitting<-CNV.fit(object.groups[i], object.controls, anno)
        raw_name<-names(fitting)
        proc_name<-gsub("[^A-Za-z0-9]", "_", raw_name)
        print(proc_name)
        binning<- CNV.bin(fitting)
        detailing<- CNV.detail(binning)

        segmenting <- CNV.segment(detailing)

        label<-toString(pheno[pheno[identifier]==names(object.groups)[i], ][labels])


       if (sample_plots) {

                  if (save)  dir.create(file.path(path,label,"samples"), showWarnings = FALSE, recursive=TRUE)
                  if (save)  png(file=file.path(path,label,"samples",paste(proc_name,".png")),width = 1920, height = 480)

               CNV.genomeplot(segmenting)

               dev.off()

                 if (save)  dir.create(file.path(path,label,"samples","tables"), showWarnings = FALSE,recursive=TRUE)
                 if (save)  CNV.write(segmenting,file=file.path(path,label,"samples","tables",paste(proc_name,"segments.txt")), what = 'segments')
                 if (save)  CNV.write(segmenting, file=file.path(path,label,"samples","tables",paste(proc_name,"detail.txt")), what = 'detail')
                 if (save)  CNV.write(segmenting, file=file.path(path,label,"samples","tables",paste(proc_name,"bins.txt")), what = 'bins')
                 if (save)  CNV.write(segmenting, file=file.path(path,label,"samples","tables",paste(proc_name,"probes.txt")), what = 'probes')

                 if (save)  dir.create(file.path(path,label,"samples","plots","detail"), recursive=TRUE, showWarnings = FALSE)
                 if (save)  png(file=file.path(path,label,"samples","plots","detail",paste(proc_name,".png")),width = 1920, height = 480)

               CNV.detailplot_wrap(segmenting)
               dev.off()

               for (p in 1:length(list_of_genes)) {
                     print(list_of_genes[[p]])
                          if (save)  dir.create(file.path(path,label,"samples","plots","detail_regions"), showWarnings = FALSE, recursive=TRUE)
                          if (save) png(file=file.path(path,label,"samples","plots","detail_regions",paste(proc_name,list_of_genes[[p]],".png")),width = 1920, height = 480)

                      CNV.detailplot(segmenting, name = list_of_genes[[p]])
                      dev.off()

               }

               chr_rle<-seqnames(detail_region)
               chr_list<-as.vector(unique(chr_rle))
               for (p in chr_list) {
                     print(p)
                          if (save)  dir.create(file.path(path,label,"samples","plots","chr_details"), showWarnings = FALSE, recursive=TRUE)
                          if (save)   png(file=file.path(path,label,"samples","plots","chr_details",paste(proc_name,p,".png")),width = 1920, height = 480)

                      CNV.genomeplot(segmenting, chr = p)
                      dev.off()

               }



       }


        cnv_list[[proc_name]]<-segmenting

        colnames_list[[i]]<-binning@name

        intensityvals_list[[i]]<-binning@bin$ratio-binning@bin$shift

        }



intensity_vals<- data.frame(matrix(unlist(intensityvals_list), nrow=length(intensityvals_list[[i]])))
rownames(intensity_vals)<- names(intensityvals_list[[i]])
colnames(intensity_vals)<- colnames_list


elename_list<-list()
intensityvalgroups_list<-list()
cnv_val_list<-list()
gl_freq_list<-list()
cnvss_list<-list()

for (r in seq_along(interest_groups)) {

      print(interest_groups[[r]])

      elenames<-list()
      intensityvalgroups<-list()
      cnv_vals<-list()
      gl_freqs<-list()
      cn_list<-list()

      el_groups <- pheno[grep(paste(interest_groups[[r]],collapse='|'), pheno[[labels]]),]

      el_keep <- match(el_groups[[identifier]],names(object.groups) )


      el_intensity_vals <- intensity_vals[,c(el_keep)]
      el_cnv_list <-cnv_list[c(el_keep)]

      x1 <- DNAcopy::CNA(genomdat = el_intensity_vals,
                      chrom = as.vector(seqnames(cnv_list[[1]]@anno@bins)), maploc = values(cnv_list[[1]]@anno@bins)$midpoint,
                      data.type = "logratio", sampleid = colnames(el_intensity_vals))

      x11 <- DNAcopy::smooth.CNA(x1)


      x2 <- DNAcopy::segment(x11, verbose = 1,alpha=0.001, undo.splits="sdundo", undo.SD=2)


      glfreq<-DNAcopy::glFrequency(x2, cutoff)
      glfreq<-glfreq[mixedorder(glfreq$chrom), ]

      cnvss_list[[r]] <- el_cnv_list
      cnv_val_list[[r]] <- x2 #cnv_vals
      gl_freq_list[[r]] <- glfreq#gl_freqs
      intensityvalgroups_list[[r]] = el_intensity_vals #intensityvalgroups
      elename_list[[r]] = interest_groups[[r]] #elenames
}


names(cnvss_list)<-interest_groups
names(cnv_val_list)<-interest_groups
names(gl_freq_list)<-interest_groups
names(intensityvalgroups_list)<-interest_groups
names(elename_list)<-interest_groups


object <- new("CNV.summaryanalysis")
object@intensity_vals <- intensityvalgroups_list
object@cnv_seg_data <- cnvss_list
object@cnvs <- cnv_val_list
object@gl_freqs <- gl_freq_list

object@date <- date()
object@names <- interest_groups

if (summary_plots) {

CNV.summaryplot(object=object,
                chr = chr,
                chrX = chrX,
                chrY = chrY,
                centromere = centromere,
                main = main,
                ylim = ylim,
                set_par = set_par,
                save=save,
                path=path)

}




return (object)

})
