#!/usr/bin/Rscript 

library(boot)
library(tibble)
library(dplyr)
library(reshape2)
library(ggplot2)
library(argparse)

subsetTypes <- function(cmat,keepers){

  orinames <- colnames(cmat)
  xnames <- colnames(cmat)

  ## patterns <- c("\\.",
  ##             ",",
  ##             "_",
  ##             " ",
  ##             "-"
  ##             )

  ## for (ii in 1:length(patterns)) {
  ##   xnames <- gsub(patterns[ii],"",xnames)
  ##   keepers <- gsub(patterns[ii],"",keepers)
  ## }

  ynames <- c()

  for (ii in 1:length(keepers)) {

    idxs <- as.vector(sapply(tolower(keepers[ii]),
                             grepl,tolower(xnames)))
    if (sum(idxs) > 0) {
      ynames <- c(ynames,orinames[idxs])
      }
   }

  ynames <- unique(ynames)


   if (length(ynames) > 0) {
      cmat <- cmat[orinames,ynames]
   } else {
     print("[ERROR] : Bad subsetting")
   }

   return(cmat)
}


bootCorr <- function(data,
                     replicates = 100,
                     alpha = 0.05
                     ) {

  prepArr <- function() {
    z <- array(0,dim = c(ncol(data),
                         ncol(data), replicates ))
    return(z)
  }

  setNames <- function(x) {
    colnames(x) <- colnames(data)
    rownames(x) <- colnames(data)
    return(x)
  }


  nsamples <- nrow(data)
  corrVals <- prepArr()

  for (ii in 1:replicates) {
    corrVals[,,ii] <- cor(data %>% sample_n(size = nsamples,
                                            replace = T),
                          method = "pearson")

  }

  res <- list()

  res$mean <- apply(corrVals,
                    c(1,2),
                    mean)

  res$ci.lower <- apply(corrVals,
                        c(1,2),
                        quantile,
                        alpha / 2.0)
  res$ci.upper <- apply(corrVals,
                        c(1,2),
                        quantile,
                        1-alpha / 2.0)
  
  res$outside <- ifelse(res$ci.lower * res$ci.upper > 0,
                        1,
                        0)

  res <- lapply(res,setNames)

  return(res)

}

parser <- ArgumentParser()

parser$add_argument("-fp","--proportion_files",
                    nargs = '+',
                    default = NULL)
parser$add_argument("-dp",
                    "--proportion_dirs",
                    nargs = '+',
                    default = NULL)

parser$add_argument("-o","--out_dir",
                    default = "/tmp")

parser$add_argument("-r",
                    "--replicates",
                    default = 1000,
                    type = 'integer')

parser$add_argument("-t",
                    "--tag",
                    default = '',
                    type = 'character')

parser$add_argument("-s",
                    "--subset",
                    nargs = '+',
                    default = NULL)

parser$add_argument("-nm",
                    "--name_map",
                    default = NULL)



parser$add_argument("-or","--order",
                    default = NULL,
                    type = 'character',
                    help = paste0('file containing ordering of types. ',
                                  'write each cell type to include on ',
                                  'a new line.')
                    )


args <- parser$parse_args()


if (!(is.null(args$proportion_dirs))) {
  dname <- args$proportion_dirs 
  pths <- list.files(dname,
                    recursive = T)

  crit1 <- pths %>%
    grepl(pattern = "tsv")

  crit2 <- pths %>%
    basename() %>%
    grepl(pattern = "^W") 

  crit <- crit1 & crit2
  pths <- file.path(dname,pths[crit])
}

if (!(is.null(args$proportion_files))) {
  pths <- args$proportion_files[file.exists(args$proportion_files)]
}

wmats <- list()
types <- c()

for (res in seq_along(pths)) {

  wmats[[res]] <- read.table(pths[res],
                             sep = '\t',
                             header = T,
                             row.names = 1,
                             check.names = F,
                             )
  ifelse(length(types) > 0,
         types <- intersect(types,
                   colnames(wmats[[res]])),
         types <- colnames(wmats[[res]])
         )

  rownames(wmats[[res]]) <- paste(res,
                                  rownames(wmats[[res]]),
                                  sep = '_')

}

wmats <- lapply(wmats,
                function(x){select(x,
                                   types)})

wmat <- do.call("rbind",
                wmats)

#print(colnames(wmat))

remove(wmats)

if (!(is.null(args$order))) {
  newOrder <- read.table(args$order,
                         sep = '\n')
  newOrder <- as.vector(newOrder[,1])
  newOrder <- gsub(" ",'-',newOrder)
  colnames(wmat) <- gsub("\\."," ",colnames(wmat))
  colnames(wmat) <- gsub("  | ","-",colnames(wmat))
  newOrder <- intersect(newOrder,colnames(wmat))

  if (length(newOrder) > 1) {
    wmat <- wmat[,newOrder]
    print("readjusting order")
  } else {
    print("none of the specified types present")
    print("no readjustment of order")
  }
}

res <- bootCorr(wmat,
                replicates = args$replicates)

res$mean <- as.matrix(res$mean)

diag(res$mean) <- NA


res$mean <- as.data.frame(res$mean)

if (!is.null(args$subset)) {
  res$mean <- subsetTypes(res$mean,
                          keepers = args$subset)

  res$outside <- res$outside[rownames(res$mean),colnames(res$mean)]
}


if (!is.null(args$name_map)) {
  reffer <- read.table(args$name_map,
                       sep = ',',
                       header=T,
                       row.names = 1,
                       check.names = F,
                       )

  ## rownames(reffer) <-  gsub(" ","\\.",rownames(reffer))
  ## rownames(reffer) <-  gsub("-","\\.",rownames(reffer))
  ## rownames(reffer) <-  gsub("\\+","\\.",rownames(reffer))


  mapper <- as.vector(reffer$new)
  names(mapper) <- rownames(reffer)


  oCn <- colnames(res$mean)
  oRn <- rownames(res$mean)


  nCn <- mapper[oCn]
  nRn <- mapper[oRn] 


  ordrC <- as.vector( mapper[mapper %in% nCn] )
  ordrR <- as.vector( mapper[mapper %in% nRn] )

  print(ordrC)

  colnames(res$mean) <- nCn
  rownames(res$mean) <- nRn

  res$mean <- res$mean[ordrR,ordrC]

  colnames(res$outside) <- nCn
  rownames(res$outside) <- nRn

  res$outside <- res$outside[ordrR,ordrC]

}

res$mean <- res$mean[,rev(colnames(res$mean))]
res$outside <- res$outside[,rev(colnames(res$outside))]

long_res <- res$mean %>%
  as.matrix() %>%
  melt()

long_res$sig <- res$outside %>%
  as.matrix %>%
  melt() %>%
  .$value %>%
  as.numeric()

colnames(long_res) <- c("type1",
                        "type2",
                        "av.corr",
                        "sig")

if (!(is.null(args$tag))) {
  oname <- paste(args$tag,"corr-plot.png",sep = '-')
} else {
  oname <- "corr-plot.png"
}

n_types <- max(dim(res$mean))

png(file.path(args$out_dir,oname),
    width = 2000 / 22 * n_types,
    height = 2000 / 22 * n_types,
    units = "px")

ggplot(long_res, aes(type1,
                     type2,
                     fill = av.corr,
                     color = sig)
       ) +

  geom_tile(aes(width = 0.9,
                height = 0.9),
            vjust = 0,
            hjust = 0,
            size = 1) +

  coord_equal() +

  scale_color_gradient(low = "gray50",
                       high = "black",
                       guide = FALSE
                       ) +

  scale_fill_gradient2(low = "#15007e",
                       mid = "white",
                       high = "#a12727",
                      space = "Lab",
                      na.value = "grey50",
                      guide = "colourbar",
                      aesthetics = "fill") + 
  xlab("") +
  ylab("") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5),
        text =element_text(size=40),
        panel.background = element_rect("white", "white"),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15)
        )

dev.off()
