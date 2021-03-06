library(stringi)
library(animalr)
library(rmsutilityr)
library(tacr)
library(nprcmanager)
library(RODBC)

rhesus_file <- "../inst/extdata/BreedingGroups1_4MendozaTest.csv"
rh_df <- read.csv(file = rhesus_file, header = TRUE, sep = ",", 
                  stringsAsFactors = FALSE)
rh_df$mamu1 <- stri_c(rh_df$Mamu.A.Haplotype.1, "_", 
                      rh_df$Mamu.B.Haplotype.1)
rh_df$mamu2 <- stri_c(rh_df$Mamu.A.Haplotype.2, "_", 
                      rh_df$Mamu.B.Haplotype.2)
new_rh_df <- rh_df[ , c(1, 7, 8, 2)]
names(new_rh_df) <- c("id", "mamu1", "mamu2", "sex")
female_breeders <- data.frame(id = new_rh_df$id[new_rh_df$sex == "F"])
female_breeders$id <- as.character(female_breeders$id)
write.csv(female_breeders, file = "../reports/2017_rhesus_female_breeders.csv",
          row.names = FALSE)

conn <- odbcConnect("frogstar-vortex-animal-sa")
ped <- get_direct_ancestors(conn, new_rh_df$id)
ped <- add_birth_date(conn, ped)
ped <- add_death_date(conn, ped)
ped <- add_departure_date(conn, ped)
names(ped)[names(ped) == "date"] <- "departure"
ped$description <- NULL
ped$dispositionType <- NULL

ped$birth_date <- format(ped$birth_date, format = "%Y-%m-%d")
ped$death_date <- format(ped$death_date, format = "%Y-%m-%d")
ped$departure <- format(ped$departure, format = "%Y-%m-%d")
ped$id <- stri_trim_both(ped$id)
ped$sire_id <- stri_trim_both(ped$sire_id)
ped$dam_id <- stri_trim_both(ped$dam_id)
ped_qc <- qcStudbook(ped)

write.csv(ped_qc, file = "../reports/2017_rhesus_breeding_groups_ped.csv",
          row.names = FALSE)
#bab_df <- read.csv(file = "/Users/msharp/Documents/Development/R/r_workspace/library/nprcmanager/inst/extdata/baboon_breeders_ped_genotype.csv", 
#                   header = TRUE, sep = ",", 
#                   stringsAsFactors = FALSE)
#bab_qc <- qcStudbook(bab_df)
#p <- bab_qc
p <- ped_qc
p["ped.num"] <- findPedigreeNumber(p$id, p$sire, p$dam)
p["gen"] <- findGeneration(p$id, p$sire, p$dam)
gv <- reportGV(p, gu.iter = 10,
               gu.thresh = 0L,
               by.id = TRUE,
               updateProgress = NULL)
## ignore can be set to list(c("F", "F")), but this causes all female 
## relatedness to be ignored which does not make sense for breeding group 
## formation.
## 
bab_breeders <- c("1X0872",
                  "1X0945",
                  "1X0951",
                  "1X1032",
                  "1X1125",
                  "1X1126",
                  "1X1146",
                  "1X1152",
                  "1X1155",
                  "1X1181",
                  "1X1237",
                  "1X1288",
                  "1X1487",
                  "1X1672",
                  "1X1765",
                  "1X1773",
                  "1X1939",
                  "1X1947",
                  "1X2055",
                  "1X2645")
numGp <- 2L
grp <- groupAssign(female_breeders$id, gv[["kinship"]], p, 
                   threshold = 0.0625, ignore = NULL, 
                   min.age = 0L, iter = 3, numGp = numGp, 
                   updateProgress = NULL, withKin = TRUE)
i <- 1
grp$group[[i]]; grp$score; grp$groupKin[[i]]
i <- 2
grp$group[[i]]; grp$score; grp$groupKin[[i]]
i <- 3
grp$group[[i]]; grp$score; grp$groupKin[[i]]
i <- 4
grp$group[[i]]; grp$score; grp$groupKin[[i]]
i <- 5
grp$group[[i]]; grp$score; grp$groupKin[[i]]
kmat <- filterKinMatrix(female_breeders$id, gv[["kinship"]])
print(kmat)
