#!/usr/bin/Rscript

library(gridExtra)
library(argparse)

prs <- ArgumentParser()

prs$add_argument("-f",
                 "--file_name",
                 help ="full path to file",
                 required = T)

prs$add_argument("-un",
                 "--use_numbers",
                 help ="add left column with number of item ",
                 default = F,
                 action = 'store_true'
                 )

prs$add_argument("-nt","--n_top",
                 help ="number of top genes to select",
                 default = NULL,
                 type = 'double')

prs$add_argument("-o",
                 "--out_name",
                 help ="full path to output",
                 default = NULL)

args <- prs$parse_args()

if (is.null(args$out_name)) {
  args$out_name <- "top-genes.png"
}

#pth <- "/tmp/tls/tls-associated.tsv"

pth <- args$file_name

data <- read.table(pth,
                  header = 1,
                  row.names = 1,
                  sep = '\t')

data$gene <- rownames(data)

if (args$use_numbers ) {
  data$rank <- c(1:nrow(data))
}

data <- data[,c(3,2,1)]
colnames(data) <- c("Rank","Gene","Coefficient")
data$Coefficient <- round(data$Coefficient,5)


png(filename = args$out_name,
    width=500,height=1000,bg = "white")

thm <- ttheme_default(base_size = 30,
                      base_family = 'calibri')

grid.table(data[1:args$n_top,],
           rows = NULL,
           cols = colnames(data),
           theme = thm)

dev.off()
