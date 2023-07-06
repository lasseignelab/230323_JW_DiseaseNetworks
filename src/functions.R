<<<<<<< HEAD
# function ppi construction
ppi <- function(species_code, score_cutoff, string_version, data){
  print("creating a new STRINGdb class with specified version, species, and score threshold")
  string_db <- STRINGdb$new(
    version = as.character(string_version),
    species = as.numeric(species_code),
    score_threshold = score_cutoff
  ) #
  # create a TF data frame with one column
  TF <- data.frame(data[, 1]) # grabs the TF column from this data
  # change the colname to "TF"
  colnames(TF) <- "TF"
  TF <- unique(TF)
  
  print("mapping the TFs to STRINGdb dataset")
  TF_mapped <- string_db$map(TF,
                             "TF",
                             removeUnmappedRows = TRUE
  ) # Warning:  we couldn't map to STRING 0% of your identifiers. Setting removeUnmappedRows to TRUE or FALSE doesn't make a difference.
  # collect the interactions between the TF of interest
  ppi_tmp <- string_db$get_interactions(TF_mapped$STRING_id)[, c(1, 2, 3)] # contains duplicate data
  ppi_tmp <- unique(ppi_tmp) # removing duplicates
  
  # store the PPI by using original identifier.
  output <- data.frame(
    from = TF_mapped[match(ppi_tmp$from, TF_mapped$STRING_id), 1],
    to = TF_mapped[match(ppi_tmp$to, TF_mapped$STRING_id), 1],
    score = ppi_tmp$combined_score
  )
  y <- nrow(output) 
  y
  
  x <- head(output) # creates a dataframe with a "from" and a "to" column for PPI
  x 
  
  return(output)
}

# function for mouse ppi construction 
mus_ppi <- function(species_code, score_cutoff, string_version, data){
  print("creating a new STRINGdb class with specified version, species, and score threshold")
  string_db <- STRINGdb$new(
    version = as.character(string_version),
    species = as.numeric(species_code),
    score_threshold = score_cutoff
  ) #
  # create a TF data frame with one column
  TF <- data.frame(data[, 1]) # grabs the TF column from this data
  # change the colname to "TF"
  colnames(TF) <- "TF"
  TF <- unique(TF)
  
  print("mapping the TFs to STRINGdb dataset")
  TF_mapped <- string_db$map(TF,
                             "TF",
                             removeUnmappedRows = TRUE
  ) # Warning:  we couldn't map to STRING 0% of your identifiers. Setting removeUnmappedRows to TRUE or FALSE doesn't make a difference.
  # collect the interactions between the TF of interest
  ppi_tmp <- string_db$get_interactions(TF_mapped$STRING_id)[, c(1, 2, 3)] # contains duplicate data
  ppi_tmp <- unique(ppi_tmp) # removing duplicates
  
  # store the PPI by using original identifier.
  output <- data.frame(
    from = TF_mapped[match(ppi_tmp$from, TF_mapped$STRING_id), 1],
    to = TF_mapped[match(ppi_tmp$to, TF_mapped$STRING_id), 1],
    score = ppi_tmp$combined_score
  )
  y <- nrow(output) 
  y
  
  x <- head(output) # creates a dataframe with a "from" and a "to" column for PPI
  x 
  
  
  # All protein names returned are in all capital, even though they are mouse (as specified by 10090 above). R is case sensitive, so need to convert these to be lowercase with first letter capital
  from <- str_to_title(output$from)
  to <- str_to_title(output$to)
  score <- str_to_title(output$score)
  
  output <- data.frame(from, to, score)
  
  return(output)
}

# Function-scatterplot
scatterPlot <- function(
    simMatrix, reducedTerms, size = "score", addLabel = TRUE,
    labelSize = 3) {
  if (!all(sapply(c("ggplot2", "ggrepel"), requireNamespace,
                  quietly = TRUE
  ))) {
    stop("Packages ggplot2, ggrepel and/or its dependencies not available. ",
         "Consider installing them before using this function.",
         call. = FALSE
    )
  }
  x <- cmdscale(as.matrix(as.dist(1 - simMatrix)),
                eig = TRUE,
                k = 2
  )
  df <- cbind(as.data.frame(x$points), reducedTerms[match(
    rownames(x$points),
    reducedTerms$go
  ), c("term", "parent", "parentTerm", "size")])
  p <- ggplot2::ggplot(df, ggplot2::aes(x = "V1", y = "V2", color = "parentTerm")) +
    ggplot2::geom_point(ggplot2::aes(size = size), alpha = 0.5) +
    ggplot2::scale_color_discrete(guide = "none") +
    ggplot2::scale_size_continuous(
      guide = "none",
      range = c(0, 25)
    ) +
    ggplot2::scale_x_continuous(name = "") +
    ggplot2::scale_y_continuous(name = "") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank()
    )
  if (addLabel) {
    p + ggrepel::geom_label_repel(aes(label = "parentTerm"),
                                  data = subset(df, parent == rownames(df)), box.padding = grid::unit(
                                    0.5,
                                    "lines"
                                  ), size = labelSize, max.overlaps = 20
    )
  } else {
    p
  }
}


# aggregate matrix function
##Matrix.utils function not loading into Docker, so grabbed source code (https://rdrr.io/cran/Matrix.utils/src/R/Matrix.utils.R)
aggregate.Matrix<-function(x,groupings=NULL,form=NULL,fun='sum',...)
{
  if(!is(x,'Matrix'))
    x<-Matrix(as.matrix(x),sparse=TRUE)
  if(fun=='count')
    x<-x!=0
  groupings2<-groupings
  if(!is(groupings2,'data.frame'))
    groupings2<-as(groupings2,'data.frame')
  groupings2<-data.frame(lapply(groupings2,as.factor))
  groupings2<-data.frame(interaction(groupings2,sep = '_'))
  colnames(groupings2)<-'A'
  if(is.null(form))
    form<-as.formula('~0+.')
  form<-as.formula(form)
  mapping<-dMcast(groupings2,form)
  colnames(mapping)<-substring(colnames(mapping),2)
  result<-t(mapping) %*% x
  if(fun=='mean')
    result@x<-result@x/(aggregate.Matrix(x,groupings2,fun='count'))@x
  attr(result,'crosswalk')<-grr::extract(groupings,match(rownames(result),groupings2$A))
  return(result)
}

# dmCast function
dMcast<-function(data,formula,fun.aggregate='sum',value.var=NULL,as.factors=FALSE,factor.nas=TRUE,drop.unused.levels=TRUE)
{
  values<-1
  if(!is.null(value.var))
    values<-data[,value.var]
  alltms<-terms(formula,data=data)
  response<-rownames(attr(alltms,'factors'))[attr(alltms,'response')]
  tm<-attr(alltms,"term.labels")
  interactionsIndex<-grep(':',tm)
  interactions<-tm[interactionsIndex]
  simple<-setdiff(tm,interactions)
  i2<-strsplit(interactions,':')
  newterms<-unlist(lapply(i2,function (x) paste("paste(",paste(x,collapse=','),",","sep='_'",")")))
  newterms<-c(simple,newterms)
  newformula<-as.formula(paste('~0+',paste(newterms,collapse='+')))
  allvars<-all.vars(alltms)
  data<-data[,c(allvars),drop=FALSE]
  if(as.factors)
    data<-data.frame(lapply(data,as.factor))
  characters<-unlist(lapply(data,is.character))
  data[,characters]<-lapply(data[,characters,drop=FALSE],as.factor)
  factors<-unlist(lapply(data,is.factor))
  #Prevents errors with 1 or fewer distinct levels
  data[,factors]<-lapply(data[,factors,drop=FALSE],function (x) 
  {
    if(factor.nas)
      if(any(is.na(x)))
      {
        levels(x)<-c(levels(x),'NA')
        x[is.na(x)]<-'NA'
      }
    if(drop.unused.levels)
      if(nlevels(x)!=length(na.omit(unique(x))))
        x<-factor(as.character(x))
    y<-contrasts(x,contrasts=FALSE,sparse=TRUE)
    attr(x,'contrasts')<-y
    return(x)
  })
  #Allows NAs to pass
  attr(data,'na.action')<-na.pass
  result<-sparse.model.matrix(newformula,data,drop.unused.levels = FALSE,row.names=FALSE)
  brokenNames<-grep('paste(',colnames(result),fixed = TRUE)
  colnames(result)[brokenNames]<-lapply(colnames(result)[brokenNames],function (x) {
    x<-gsub('paste(',replacement='',x=x,fixed = TRUE) 
    x<-gsub(pattern=', ',replacement='_',x=x,fixed=TRUE) 
    x<-gsub(pattern='_sep = \"_\")',replacement='',x=x,fixed=TRUE)
    return(x)
  })
  
  result<-result*values
  if(isTRUE(response>0))
  {
    responses=all.vars(terms(as.formula(paste(response,'~0'))))
    result<-aggregate.Matrix(result,data[,responses,drop=FALSE],fun=fun.aggregate)
  }
  return(result)
}

# Functions for PANDA network job for loop
#make function for loading .Rdata to reassign to variable in loop 
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}


#run 'makePanda' function to make PANDA networks
makePanda <- function(motif, ppi, expression){
  
  #remove expression values = 0 
  #expression <- expression[rowSums(expression[])>0,] #removing values with zeroes, this is commented out because removing them here will cause different sized matrices and result in downstream network comparisons not being possible
  
  #run PANDA
  panda(expr = expression, ppi = ppi, motif = motif, progress = TRUE, mode = "intersection") #running in default mode 
  
}

# Functions for Setbp1_AllCortex_PANDAComparison_positive_01

## function-for-functional-enrichment-pathway-analysis
fea <- function(genes, organism, max_term = 1000, min_term = 5){
  # create gprofiler2 query ---
  fea_result <- gost(query = genes, organism = organism, ordered_query = FALSE, multi_query = FALSE, significant = 
                       TRUE, exclude_iea = FALSE, measure_underrepresentation = FALSE, evcodes = TRUE, user_threshold = 
                       0.05, correction_method = "bonferroni", domain_scope = "annotated", numeric_ns =
                       "", sources = NULL, as_short_link = FALSE) 
  # remove arbitrary pathways --- do not want pathways too "generic"
  fea_result_filt <- fea_result$result %>% dplyr::filter(., term_size < max_term & term_size > min_term) #std for max_term is 1000 and min_term is 10
  # # select the top ___ pathways for plotting --- suggest 50, could do more, but difficult to see visually
  # fea_result_filt <- fea_result %>% top_n(n = pathway_number)
  return(fea_result_filt)
}

fea_no_sig <- function(genes, organism, max_term = 1000, min_term = 5){
  # create gprofiler2 query ---
  fea_result <- gost(query = genes, organism = organism, ordered_query = FALSE, multi_query = FALSE, significant = 
                       FALSE, exclude_iea = FALSE, measure_underrepresentation = FALSE, evcodes = TRUE, user_threshold = 
                       0.05, correction_method = "bonferroni", domain_scope = "annotated", numeric_ns =
                       "", sources = NULL, as_short_link = FALSE) 
  #fea_result_filt <- fea_result
  # remove arbitrary pathways --- do not want pathways too "generic"
  #\fea_result_filt <- fea_result$result %>% dplyr::filter(., term_size < max_term & term_size > min_term) #std for max_term is 1000 and min_term is 10
  # # select the top ___ pathways for plotting --- suggest 50, could do more, but difficult to see visually
  # fea_result_filt <- fea_result %>% top_n(n = pathway_number)
  return(fea_result_filt)
}

## function-for-plotting-functional-enrichment-pathway-analysis
bubbleplot <- function(fea_result_filt){
  plot <- ggplot(fea_result_filt, aes(x = intersection_size, y = reorder(term_name, -p_value), size = recall, fill =
                                        p_value)) +
    geom_point(alpha = 0.7, shape = 21) +
    scale_size(range = c(2, 10), name = "# Genes Matched to Term") + 
    scale_fill_distiller(palette = "Purples") + 
    labs(x = "Intersection Size", y = "Functional Enrichment Terms") +
    theme_minimal()
  return(plot)
}


## function-GO-positive-negative-combined
combined_fea <- function(genes, organism, max_term = 1000, min_term = 5){
  posgenes <- as.list(genes %>% filter(value == "positive"))
  pos_fea_filt <- fea(genes= posgenes, organism = organism, max_term = max_term, min_term = min_term) %>% mutate(direction = "positive")
  neggenes <- as.list(genes %>% filter(value == "negative"))
  neg_fea_filt <- fea(genes = neggenes, organism = organism, max_term = max_term, min_term = min_term) %>% mutate(direction = "negative")
  combined_filt <- rbind(pos_fea_filt, neg_fea_filt)
  combined_fea_result <- list(control_enriched = neg_fea_filt, condition_enriched = pos_fea_filt, combined = combined_filt)
  return(combined_fea_result)
}

#function DGE fea
fea_DGE <- function(genes, organism){
  # create gprofiler2 query ---
  fea_result <- gost(query = genes, organism = organism, ordered_query = FALSE, multi_query = FALSE, significant = 
                       TRUE, exclude_iea = FALSE, measure_underrepresentation = FALSE, evcodes = TRUE, user_threshold = 
                       0.05, correction_method = "bonferroni", domain_scope = "annotated", numeric_ns =
                       "", sources = NULL, as_short_link = FALSE) 
  # remove arbitrary pathways --- do not want pathways too "generic"
  fea_result <- fea_result$result %>% filter(term_size < 1000 | term_size > 10)
  # select the top 50 pathways for plotting --- could do more, but difficult to see visually
  fea_result_filt <- fea_result %>% top_n(n = 50)
  return(fea_result_filt)
}

combined_fea_DGE <- function(genes, organism){
  upgenes <- as.list(genes %>% filter(diffexpressed == "UP"))
  up_fea_filt <- fea(genes = upgenes, organism = organism) %>% mutate(direction = "upregulated")
  downgenes <- as.list(genes %>% filter(diffexpressed == "DOWN"))
  down_fea_filt <- fea(genes = downgenes, organism = organism) %>% mutate(direction = "downregulated")
  combined_filt <- rbind(up_fea_filt, down_fea_filt)
  fea_result <- list(downregulated = down_fea_filt, upregulated = up_fea_filt, combined = combined_filt)
  return(fea_result)
}

#setbp1 gene set enrichment functions side by side fea
fea_set <- function(genes, organism){
  # create gprofiler2 query ---
  fea_result <- gost(query = genes, organism = organism, ordered_query = FALSE, multi_query = FALSE, significant = 
                       FALSE, exclude_iea = FALSE, measure_underrepresentation = FALSE, evcodes = TRUE, user_threshold = 
                       0.05, correction_method = "bonferroni", domain_scope = "annotated", numeric_ns =
                       "", sources = NULL, as_short_link = FALSE) 
  # remove arbitrary pathways --- do not want pathways too "generic"
  #fea_result <- fea_result$result %>% filter(term_size < 1000 | term_size > 10) $removed because gene set is small
  # select the top 50 pathways for plotting --- could do more, but difficult to see visually
  fea_result_filt <- fea_result$result %>% top_n(n = 50)
  return(fea_result_filt)
}

combined_fea_set <- function(genes, organism){
  upgenes <- as.list(genes %>% filter(diffexpressed == "UP"))
  up_fea_filt <- fea_set(genes = upgenes, organism = organism) %>% mutate(direction = "upregulated")
  downgenes <- as.list(genes %>% filter(diffexpressed == "DOWN"))
  down_fea_filt <- fea_set(genes = downgenes, organism = organism) %>% mutate(direction = "downregulated")
  combined_filt <- rbind(up_fea_filt, down_fea_filt)
  fea_result <- list(downregulated = down_fea_filt, upregulated = up_fea_filt, combined = combined_filt)
  return(fea_result)
}


## function-aggregate-matrix
aggregate.Matrix<-function(x,groupings=NULL,form=NULL,fun='sum',...)
{
  if(!is(x,'Matrix'))
    x<-Matrix(as.matrix(x),sparse=TRUE)
  if(fun=='count')
    x<-x!=0
  groupings2<-groupings
  if(!is(groupings2,'data.frame'))
    groupings2<-as(groupings2,'data.frame')
  groupings2<-data.frame(lapply(groupings2,as.factor))
  groupings2<-data.frame(interaction(groupings2,sep = '_'))
  colnames(groupings2)<-'A'
  if(is.null(form))
    form<-as.formula('~0+.')
  form<-as.formula(form)
  mapping<-dMcast(groupings2,form)
  colnames(mapping)<-substring(colnames(mapping),2)
  result<-t(mapping) %*% x
  if(fun=='mean')
    result@x<-result@x/(aggregate.Matrix(x,groupings2,fun='count'))@x
  attr(result,'crosswalk')<-grr::extract(groupings,match(rownames(result),groupings2$A))
  return(result)
}

## DmCast Function 
dMcast<-function(data,formula,fun.aggregate='sum',value.var=NULL,as.factors=FALSE,factor.nas=TRUE,drop.unused.levels=TRUE)
{
  values<-1
  if(!is.null(value.var))
    values<-data[,value.var]
  alltms<-terms(formula,data=data)
  response<-rownames(attr(alltms,'factors'))[attr(alltms,'response')]
  tm<-attr(alltms,"term.labels")
  interactionsIndex<-grep(':',tm)
  interactions<-tm[interactionsIndex]
  simple<-setdiff(tm,interactions)
  i2<-strsplit(interactions,':')
  newterms<-unlist(lapply(i2,function (x) paste("paste(",paste(x,collapse=','),",","sep='_'",")")))
  newterms<-c(simple,newterms)
  newformula<-as.formula(paste('~0+',paste(newterms,collapse='+')))
  allvars<-all.vars(alltms)
  data<-data[,c(allvars),drop=FALSE]
  if(as.factors)
    data<-data.frame(lapply(data,as.factor))
  characters<-unlist(lapply(data,is.character))
  data[,characters]<-lapply(data[,characters,drop=FALSE],as.factor)
  factors<-unlist(lapply(data,is.factor))
  #Prevents errors with 1 or fewer distinct levels
  data[,factors]<-lapply(data[,factors,drop=FALSE],function (x) 
  {
    if(factor.nas)
      if(any(is.na(x)))
      {
        levels(x)<-c(levels(x),'NA')
        x[is.na(x)]<-'NA'
      }
    if(drop.unused.levels)
      if(nlevels(x)!=length(na.omit(unique(x))))
        x<-factor(as.character(x))
    y<-contrasts(x,contrasts=FALSE,sparse=TRUE)
    attr(x,'contrasts')<-y
    return(x)
  })
  #Allows NAs to pass
  attr(data,'na.action')<-na.pass
  result<-sparse.model.matrix(newformula,data,drop.unused.levels = FALSE,row.names=FALSE)
  brokenNames<-grep('paste(',colnames(result),fixed = TRUE)
  colnames(result)[brokenNames]<-lapply(colnames(result)[brokenNames],function (x) {
    x<-gsub('paste(',replacement='',x=x,fixed = TRUE) 
    x<-gsub(pattern=', ',replacement='_',x=x,fixed=TRUE) 
    x<-gsub(pattern='_sep = \"_\")',replacement='',x=x,fixed=TRUE)
    return(x)
  })
  
  result<-result*values
  if(isTRUE(response>0))
  {
    responses=all.vars(terms(as.formula(paste(response,'~0'))))
    result<-aggregate.Matrix(result,data[,responses,drop=FALSE],fun=fun.aggregate)
  }
  return(result)
}

# split violin plot function
library(ggplot2)
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin,
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                             
                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                         1))
                               quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                               aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                             } else {
                               ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                             }
                           }
)

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ...,
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE,
                              show.legend = NA, inherit.aes = TRUE) {
  layer(
    data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...)
  )
}

# function targeting-Calc on panda regNet for gene and TF
targetingCalc <- function(regNetmatrix, variable_name, edge_weight_name, condition){
  #rearrange dataframe 
  regNetmatrix <- reshape2::melt(regNetmatrix, varnames = c("TF", "gene"), value.name = "edge_weight_name")#melting dataframe
  print("dataframe melted")
  regNetmatrix$edge_weight_name_pos <- ifelse(regNetmatrix$edge_weight_name < 0, 0, regNetmatrix$edge_weight_name) #replacing all negatives as a 0 and storing in new column
  print("subsetting only positive edge weights")
  regNetmatrix <- regNetmatrix[,c(1,2,4)]
  regNetmatrix
  
  #calculate gene targeting
  print("calculating gene targeting")
  Gene.targeting <- aggregate(.~gene, regNetmatrix[-1], sum) #removing TF column and calculating targeting for all edge weights and when edge weight is only positive
  #set column names based on condition
  print("renaming columns")
  if(condition == "het"){
    colnames(Gene.targeting) <- c("gene", "het_edge_weight_pos")
    print("plotting gene targeting score distribution")
    png(file = paste0(here("results/diff_targeting/"), variable_name, condition, "_GeneTargetingScoresDist.png"),
        width = 1000,
        height = 1000)
      hist(Gene.targeting$het_edge_weight_pos)
      dev.off()
  } else {
    colnames(Gene.targeting) <- c("gene", "ctrl_edge_weight_pos")
    print("plotting gene targeting score distribution")
    png(file = paste0(here("results/diff_targeting/"), variable_name, condition, "_GeneTargetingScoresDist.png"),
        width = 1000,
        height = 1000)
    hist(Gene.targeting$ctrl_edge_weight_pos)
    dev.off()
  }
  #reassign variable 
  print("assigning variable name to object")
  variable_name <- as.character(variable_name)
  assign(paste0(variable_name, "_gene_targeting_", condition), Gene.targeting, envir = .GlobalEnv)
  print("gene targeting calculation complete")
  
  #calculate TF targeting
  print("calculating TF targeting")
  TF.targeting <- aggregate(.~TF, regNetmatrix[-2], sum) #same as above but for TF instead of gene
  #set column names based on condition
  print("renaming columns")
  if(condition == "het"){
    colnames(TF.targeting) <- c("TF", "het_edge_weight_pos")
    print("plotting TF targeting score distribution")
    png(file = paste0(here("results/diff_targeting/"), variable_name, condition, "_TFTargetingScoresDist.png"),
        width = 1000,
        height = 1000)
    hist(TF.targeting$het_edge_weight_pos)
    dev.off()
  } else {
    colnames(TF.targeting) <- c("TF", "ctrl_edge_weight_pos")
    print("plotting TF targeting score distribution")
    png(file = paste0(here("results/diff_targeting/"), variable_name, condition, "_TFTargetingScoresDist.png"),
        width = 1000,
        height = 1000)
    hist(TF.targeting$ctrl_edge_weight_pos)
    dev.off()
  }
  #reassign variable
  print("assigning variable name to object")
  variable_name <- as.character(variable_name)
  assign(paste0(variable_name, "_TF_targeting_", condition), TF.targeting, envir = .GlobalEnv)
  print("TF targeting calculation complete")
}

# function-targeting_heatmap; used in targeting
targeting_heatmap <- function(annotation_colors, data, meta_colname, plot_path, rowtitle, plot_title, show_names = FALSE, width = 1000, height = 1000){
  #plotting all 
  ##grabbing metadata and annotations
  meta <- as.data.frame(colnames(data))
  colnames(meta) <- meta_colname
  rownames(meta) <- meta[,1]
  
  ##set heatmap annotations
  heat.anno = HeatmapAnnotation(df = meta, show_annotation_name = TRUE, col = annotation_colors)
  
  ##ensure column order matches annotation table
  rows <- rownames(data)
  data <- data[,rownames(meta), drop = FALSE]
  rownames(data) <- rows
  ##convert data to matrix
  mat <- as.matrix(data)
  
  ##plot heatmap 
  png(filename = plot_path,
      width = width,
      height = height)
  print(Heatmap(mat,
                col = colorRampPalette(brewer.pal(8,"Blues")) (25),
                heatmap_legend_param = list(title = "targeting score"),
                cluster_rows = TRUE,
                cluster_columns = TRUE,
                column_order = NULL,
                show_row_dend = TRUE,
                show_column_dend = TRUE,
                show_row_names = show_names,
                show_column_names = FALSE,
                use_raster = TRUE,
                raster_device = c("png"),
                bottom_annotation = NULL,
                top_annotation = heat.anno,
                column_title = plot_title, row_title = rowtitle, row_title_side = "right"))
  dev.off()
}

#functions from CosIA for gene conversion
updateGEx <- function (expr_mat, species) {
  # get rownames from expression data
  species_rowid_vec <- rownames(expr_mat)
  species_rowid_df <- as.data.frame(rownames(expr_mat)) |>
    dplyr::rename(X1 = `rownames(expr_mat)`)

  # use CoSIA to convert Bgee Ensembl IDs to species gene symbols
  converted_id <- CoSIA_ensembltosymbol(species, species_rowid_df)
  print("Check first converted ID")
  print(head(converted_id, 1))

  # match names between Ensembl id rownames and converted gene symbols
  matched_names <- converted_id[,2][match(species_rowid_vec, converted_id[,1])]
  print("Check first matched ID")
  print(head(matched_names, 1))

  # recursively replace Ensembl ids with gene symbol,
  ##unless a gene symbol is not there
  for (i in 1:length(matched_names)){
    if (! is.na(matched_names[i])) {
      rownames(expr_mat)[i] <- matched_names[i]
    } else {
      next
    }
  }

  # log2+1 transformation of tpm values
  transGEx <- log2(expr_mat + 1)

  # Getting unique
  unqGEx <- transGEx[!duplicated(rownames(transGEx)),]

  return(unqGEx)
}

CoSIA_ensembltosymbol <- function(species, gene_list){
  spec <- switch(
    species,
    "Homo_sapiens" = "h_sapiens",
    "Mus_musculus" = "m_musculus",
    "Rattus_norvegicus" = "r_norvegicus")

  tmp_cosia <- CoSIA::CoSIAn(
    gene_set = gene_list$X1,
    i_species = spec,
    o_species = c(spec),
    input_id = "Ensembl_id",
    output_ids = "Symbol",
    map_species = c(spec),
    map_tissues = c("brain","kidney"),
    mapping_tool = "annotationDBI",
    ortholog_database = "HomoloGene",
    metric_type = "DS_Gene"
  )

  tmp_cosia <- CoSIA::getConversion(tmp_cosia)
  print("returning CoSIA object...")
  return(tmp_cosia@converted_id)
}
=======
# Funtions for 01_EdgeSpecificity.Rmd

calculate_proportions <- function (gene_list,
                          data_table,
                          tissue_types) {
  # Subset data table for genes of interest
  gene_subsets <- list()
  for (i in 1:length(gene_list)) {
    gene_subset <- data_table[
      data_table$TargetGene == gene_list[[i]], 
    ]
    gene <- gene_list[[i]]
    gene_subsets[[gene]] <- gene_subset
  }
  gene_edge_counts <- list()
  # Calculate total number of edges for each gene of interest
  for (i in 1:length(names(gene_subsets))) {
    gene <- names(gene_subsets)[[i]]
    gene_edge_counts[[gene]] <- nrow(gene_subsets[[i]])
  }
  # Create a list of dataframes for each gene that is appropriately sized
  gene_tissue_dfs <- list()
  for (i in 1:length(names(gene_subsets))) {
    gene <- names(gene_subsets)[[i]]
    gene_tissue_dfs[[gene]] <- data.frame(
      matrix(
        ncol = 38, 
        nrow = gene_edge_counts[[gene]]
      )
    )
  }
  # Format the dataframes with the proper row and col names
  for (i in 1:length(names(gene_tissue_dfs))) {
    gene <- names(gene_tissue_dfs)[[i]]
    rownames(gene_tissue_dfs[[gene]]) <- gene_subsets[[gene]]$TF
    colnames(gene_tissue_dfs[[gene]]) <- tissue_types
  }
  # Assigne 1 for each edge that is specific per tissue
  for (i in 1:length(names(gene_tissue_dfs))) {
    gene <- names(gene_tissue_dfs)[[i]]
    for (j in 1:nrow(gene_tissue_dfs[[gene]])) {
      present <- unlist(strsplit(gene_subsets[[gene]]$Tissues[j], ","))
      for (k in 1:length(present)) {
        tissueForEdge <- present[k]
        gene_tissue_dfs[[gene]][j,which(colnames(gene_tissue_dfs[[gene]])==tissueForEdge)]<-1
      }
    }
  }
  
  # Assign 0 for TF-pkd gene pairs that are not present in each tissue
  for (i in 1:length(names(gene_tissue_dfs))) {
    gene_dataframe <- gene_tissue_dfs[[i]]
    gene_dataframe[is.na(gene_dataframe)] <- 0
    gene_tissue_dfs[[i]] <- gene_dataframe
  }
  # Sum the number of specific edges for each gene in each tissue
  edges_gene_tissues <- list()
  for (i in 1:length(names(gene_tissue_dfs))) {
    gene_dataframe <- gene_tissue_dfs[[i]]
    edges_gene_tissues[[i]] <- colSums(gene_dataframe)
  }
  # Calculate proportion of specific edges vs total edges per tissue
  proportions_per_gene <- list()
  for (i in 1:length(names(gene_tissue_dfs))) {
    gene <- names(gene_tissue_dfs)[[i]]
    gene_dataframe <- gene_tissue_dfs[[i]]
    gene_proportion <- edges_gene_tissues[[i]]/gene_edge_counts[[i]]
    proportions_per_gene[[gene]] <- gene_proportion
  }
  return(proportions_per_gene)
}


# Plotting proportions
plot_proportions <- function (proportions_df,
                              gene_list) {
  for (i in 1:length(names(formatted_proportions))) {
    gene <- names(formatted_proportions)[[i]]
    plot <- ggplot(
      formatted_proportions[[gene]],
      aes(
        x = Proportions,
        y = Tissues,
      )
    ) +
      geom_point(
        col = "tomato2",
        size = 3
      ) +
      xlim(0, 1) +
      labs(
        title = paste0(
          "Proportions of Specific Edges Per Tissue - ",
          gene_list[[i]]),
      )
    print(plot)
  }
}

permutation_test <- function(data_table,
                             tissue_types,
                             gene_oi,
                             nperm) {
  counter <- 0
  results <- data.frame(
    matrix(
      ncol = 1, 
      nrow = length(tissue_types)
    )
  )
  
  rownames(results) <- tissue_types
  colnames(results) <- gene_oi
  unique_targets <- unique(STable1_GTEx_TSInfo_Edges$TargetGene)
  gene_list <- c(
    gene_oi, 
    sample(
      unique_targets, 
      size = 999, 
      replace = FALSE
    )
  )
  print(length(gene_list))
  # Subset the dataframe based on the sampled gene list.
  permutation <- calculate_proportions(
    gene_list = gene_list,
    data_table = data_table,
    tissue_types = tissue_types
  )
  for (i in 1:length(tissue_types)) {
    for (j in 2:length(gene_list)){
      if (permutation[[1]][[tissue_types[i]]] <= permutation[[j]][[tissue_types[i]]]) {
        counter <- counter + 1
      }
    }
    p_val <- ((counter+1)/(nperm+1))
    results[i, 1] <- p_val
    counter <- 0
  }
  return(results)
}




# 
# formatted_proportions <- list()
# for (i in 1:length(names(proportions_df))) {
#   gene <- names(proportions_df)[[i]]
#   gene_dataframe <- as.data.frame(proportions_df[[gene]])
#   colnames(gene_dataframe) <- "Proportions"
#   gene_dataframe$Tissues <- rownames(gene_dataframe)
#   formatted_proportions[[gene]] <- gene_dataframe
# }

# Funtions for 01_EdgeSpecificity.Rmd

calculate_proportions <- function (gene_list,
                          data_table,
                          tissue_types) {
  # Subset data table for genes of interest
  gene_subsets <- list()
  for (i in 1:length(gene_list)) {
    gene_subset <- data_table[
      data_table$TargetGene == gene_list[[i]], 
    ]
    gene <- gene_list[[i]]
    gene_subsets[[gene]] <- gene_subset
  }
  gene_edge_counts <- list()
  # Calculate total number of edges for each gene of interest
  for (i in 1:length(names(gene_subsets))) {
    gene <- names(gene_subsets)[[i]]
    gene_edge_counts[[gene]] <- nrow(gene_subsets[[i]])
  }
  # Create a list of dataframes for each gene that is appropriately sized
  gene_tissue_dfs <- list()
  for (i in 1:length(names(gene_subsets))) {
    gene <- names(gene_subsets)[[i]]
    gene_tissue_dfs[[gene]] <- data.frame(
      matrix(
        ncol = 38, 
        nrow = gene_edge_counts[[gene]]
      )
    )
  }
  # Format the dataframes with the proper row and col names
  for (i in 1:length(names(gene_tissue_dfs))) {
    gene <- names(gene_tissue_dfs)[[i]]
    rownames(gene_tissue_dfs[[gene]]) <- gene_subsets[[gene]]$TF
    colnames(gene_tissue_dfs[[gene]]) <- tissue_types
  }
  # Assigne 1 for each edge that is specific per tissue
  for (i in 1:length(names(gene_tissue_dfs))) {
    gene <- names(gene_tissue_dfs)[[i]]
    for (j in 1:nrow(gene_tissue_dfs[[gene]])) {
      present <- unlist(strsplit(gene_subsets[[gene]]$Tissues[j], ","))
      for (k in 1:length(present)) {
        tissueForEdge <- present[k]
        gene_tissue_dfs[[gene]][j,which(colnames(gene_tissue_dfs[[gene]])==tissueForEdge)]<-1
      }
    }
  }
  
  # Assign 0 for TF-pkd gene pairs that are not present in each tissue
  for (i in 1:length(names(gene_tissue_dfs))) {
    gene_dataframe <- gene_tissue_dfs[[i]]
    gene_dataframe[is.na(gene_dataframe)] <- 0
    gene_tissue_dfs[[i]] <- gene_dataframe
  }
  # Sum the number of specific edges for each gene in each tissue
  edges_gene_tissues <- list()
  for (i in 1:length(names(gene_tissue_dfs))) {
    gene_dataframe <- gene_tissue_dfs[[i]]
    edges_gene_tissues[[i]] <- colSums(gene_dataframe)
  }
  # Calculate proportion of specific edges vs total edges per tissue
  proportions_per_gene <- list()
  for (i in 1:length(names(gene_tissue_dfs))) {
    gene <- names(gene_tissue_dfs)[[i]]
    gene_dataframe <- gene_tissue_dfs[[i]]
    gene_proportion <- edges_gene_tissues[[i]]/gene_edge_counts[[i]]
    proportions_per_gene[[gene]] <- gene_proportion
  }
  return(proportions_per_gene)
}


# Plotting proportions
plot_proportions <- function (proportions_df,
                              gene_list) {
  for (i in 1:length(names(formatted_proportions))) {
    gene <- names(formatted_proportions)[[i]]
    plot <- ggplot(
      formatted_proportions[[gene]],
      aes(
        x = Proportions,
        y = Tissues,
      )
    ) +
      geom_point(
        col = "tomato2",
        size = 3
      ) +
      xlim(0, 1) +
      labs(
        title = paste0(
          "Proportions of Specific Edges Per Tissue - ",
          gene_list[[i]]),
      )
    print(plot)
  }
}

permutation_test <- function(data_table,
                             tissue_types,
                             gene_oi,
                             nperm) {
  counter <- 0
  results <- data.frame(
    matrix(
      ncol = 1, 
      nrow = length(tissue_types)
    )
  )
  
  rownames(results) <- tissue_types
  colnames(results) <- gene_oi
  unique_targets <- unique(STable1_GTEx_TSInfo_Edges$TargetGene)
  gene_list <- c(
    gene_oi, 
    sample(
      unique_targets, 
      size = 999, 
      replace = FALSE
    )
  )
  print(length(gene_list))
  # Subset the dataframe based on the sampled gene list.
  permutation <- calculate_proportions(
    gene_list = gene_list,
    data_table = data_table,
    tissue_types = tissue_types
  )
  for (i in 1:length(tissue_types)) {
    for (j in 2:length(gene_list)){
      if (permutation[[1]][[tissue_types[i]]] <= permutation[[j]][[tissue_types[i]]]) {
        counter <- counter + 1
      }
    }
    p_val <- ((counter+1)/(nperm+1))
    results[i, 1] <- p_val
    counter <- 0
  }
  return(results)
}




# 
# formatted_proportions <- list()
# for (i in 1:length(names(proportions_df))) {
#   gene <- names(proportions_df)[[i]]
#   gene_dataframe <- as.data.frame(proportions_df[[gene]])
#   colnames(gene_dataframe) <- "Proportions"
#   gene_dataframe$Tissues <- rownames(gene_dataframe)
#   formatted_proportions[[gene]] <- gene_dataframe
# }
>>>>>>> 75c69bdb706354ebdf927f94217e9771cc1ca0da
