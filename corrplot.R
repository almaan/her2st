#!/usr/bin/Rscript 

library(boot)
library(tibble)
library(dplyr)
library(reshape2)
library(ggplot2)
library(argparse)


bootCorr <- function(data,
                     replicates = 100) {

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
                        0.05)
  res$ci.upper <- apply(corrVals,
                        c(1,2),
                        quantile,
                        0.95)
  
  res$outside <- ifelse(res$ci.lower * res$ci.upper > 0,
                        1,
                        0)

  res <- lapply(res,setNames)

  return(res)

}

parser <- ArgumentParser()

parser$add_argument("-fp","--proportion_files",
                    nargs = '+',
                    default = NULL,
                    help = "paths to proportion files to include")

parser$add_argument("-dp",
                    "--proportion_dirs",
                    nargs = '+',
                    default = NULL,
                    help = paste0('path to directory of proportion files. ',
                                  'All proportion files in dir will be used')
                    )

parser$add_argument("-od","--out_dir",
                    default = "/tmp",
                    help = 'output directory')

parser$add_argument("-r",
                    "--replicates",
                    default = 1000,
                    type = 'integer',
                    help = 'number of bootstrap replicates')

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

wmat <- do.call("rbind",wmats)
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


png(file.path(args$out_dir,"corr-plot.png"),
    width = 1000,
    height = 1000,
    units = "px")

ggplot(long_res, aes(type1,
                     type2,
                     fill = av.corr,
                     color = sig)
       ) +

  geom_tile(aes(width = 0.9,
                height = 0.9),
            size = 1) +

  scale_color_gradient(low = "gray50",
                       high = "black",
                       guide = FALSE
                       ) +

  scale_fill_gradient2(low = "blue",
                       mid = "white",
                      high = "red",
                      space = "Lab",
                      na.value = "grey50",
                      guide = "colourbar",
                      aesthetics = "fill") + 
  xlab("") +
  ylab("") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

dev.off()
