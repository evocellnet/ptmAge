
library(biomaRt)
library(RMySQL)

args <- commandArgs(TRUE)

comparaDb <- args[1]
outputFile <- args[2]

marts <- data.frame(hosts=c("ensembl.org",
                            "fungi.ensembl.org",
                            "metazoa.ensembl.org",
                            "protists.ensembl.org",
                            "plants.ensembl.org"),
                    marts=c("ENSEMBL_MART_ENSEMBL",
                            "fungal_mart",
                            "metazoa_mart",
                            "protist_mart",
                            "plants_mart"))

alldatasets <- apply(marts, 1, function(x) listDatasets(useMart(host=x[1],x[2])))

result <- do.call("rbind", lapply(1:length(alldatasets), function(x)
  cbind(alldatasets[[x]],
        rep(marts[x,1],nrow(alldatasets[[x]])),
        rep(marts[x,2],nrow(alldatasets[[x]])))))

names(result)[4:5] <- c("host","mart")
result$organism <- sapply(strsplit(result$description, " genes "), function(x) x[1])

#Getting taxids and merging
con <- dbConnect(MySQL(),
                 user="anonymous",
                 password="",
                 port=3306,
                 dbname=comparaDb,
                 host="ensembldb.ensembl.org")

taxidQuery <- "SELECT taxon_id,name FROM ncbi_taxa_name WHERE name_class = 'scientific name'";
rs <- dbSendQuery(con, taxidQuery)
taxidDf <- fetch(rs, n=-1)
taxidDf$name <- gsub(" Group$","",taxidDf$name)

completeresult <- merge(result, taxidDf, all.x=TRUE, by.x="organism", by.y="name")

## For the missing cases, map them using the assembly version
taxidFromGenome <- "SELECT assembly,taxon_id FROM genome_db";
rs <- dbSendQuery(con, taxidFromGenome)
taxidFromGenomeAssembly <- fetch(rs, n=-1)
assembly2taxid <- taxidFromGenomeAssembly$taxon_id
names(assembly2taxid) <- taxidFromGenomeAssembly$assembly
completeresult[is.na(completeresult$taxon_id), "taxon_id"] <- assembly2taxid[completeresult[is.na(completeresult$taxon_id), "version"]]

write.csv(completeresult,
          quote=FALSE,
          row.names=FALSE,
          file=outputFile)
