
library(SnapATAC)
library(pagoda2)
library(config)

# directory of link table
print(getwd())
config <- config::get(file = "conf/snareR.config.yaml")
link.dir <- config$link.dir
link.table <- read.table(file.path(link.dir, "links.txt"), sep="\t", header=TRUE)
saveRDS(link.table, "link.table.rds")
