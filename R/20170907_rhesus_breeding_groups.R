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
write.csv(female_breeders, file = "../reports/2017_rhesus_female_breeders.csv",
          row.names = FALSE)

conn <- odbcConnect("frogstar-vortex-animal-sa")
ped <- get_direct_ancestors(conn, new_rh_df$id)
ped <- add_birth_date(conn, ped)
ped <- add_death_date(conn, ped)
ped$birth_date <- format(ped$birth_date, format = "%Y-%m-%d")
ped$death_date <- format(ped$death_date, format = "%Y-%m-%d")
ped$id <- stri_trim_both(ped$id)
ped$sire_id <- stri_trim_both(ped$sire_id)
ped$dam_id <- stri_trim_both(ped$dam_id)
ped_qc <- qc.Studbook(ped)

write.csv(ped_qc, file = "../reports/2017_rhesus_breeding_groups_ped.csv",
          row.names = FALSE)