#!/bin/R
# reformat example inputs 

library(data.table)
library(readxl)
library(writexl)

# reformat read manifests 
adu830 = fread("~/Documents/motrpac/CLINICAL/metadata/manifests/ADU830-10060.csv", sep=",")
stanford_adu830 = read_input("~/Documents/motrpac/CLINICAL/metadata/shipment/OLD/Stanford_ADU830-10060_120720.xlsx")
ped830 = fread("~/Documents/motrpac/CLINICAL/metadata/manifests/PED830-10062.csv", sep=",")
stanford_ped830 = read_input("~/Documents/motrpac/CLINICAL/metadata/shipment/OLD/Stanford_PED830-10062_120720.xlsx")

manifest = rbindlist(list(adu830, ped830), fill=T)
shipment = rbindlist(list(stanford_adu830, stanford_ped830), fill=T)

# remove unused column names 
manifest = manifest[,.(ManifestID, Study, codedSiteID, PID, VialLabel, barcode, visitCode, SampleTypeCode, randomGroupCode, sex_psca, calculatedAge)]
shipment = shipment[,.(Box, Position, Viallabel, `2D Barcode`, `ppt type`, `Sample Type`)]

# replace . with NA
shipment[shipment=='.'] = NA
manifest[manifest=='.'] = NA

## change PID, viallabel, barcode
# PID 
pids = na.omit(unique(manifest[,PID]))
new_pids = sample(1e7:2e7, size=length(pids), replace=FALSE)
names(new_pids) = pids
manifest[!is.na(PID) ,PID := new_pids[as.character(PID)]]

# viallabel 
vl = na.omit(c(unique(manifest[,VialLabel]), unique(shipment[,Viallabel])))
new_vl = as.numeric(sample(1e10:2e10, size=length(vl), replace=FALSE))
names(new_vl) = vl
manifest[!is.na(VialLabel), VialLabel := new_vl[as.character(VialLabel)]]
shipment[!is.na(Viallabel), Viallabel := new_vl[as.character(Viallabel)]]

# barcode 
bc = na.omit(c(unique(manifest[,barcode]), unique(shipment[,`2D Barcode`])))
new_bc = as.numeric(sample(1e9:2e9, size=length(bc), replace=FALSE))
names(new_bc) = bc
manifest[!is.na(barcode), barcode := new_bc[as.character(barcode)]]
shipment[!is.na(`2D Barcode`), `2D Barcode` := new_bc[as.character(`2D Barcode`)]]

# write out 
write.table(manifest, "inst/extdata/biospecimen-metadata.csv", sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE)
write_xlsx(shipment, "inst/extdata/shipment-manifest.xlsx")
