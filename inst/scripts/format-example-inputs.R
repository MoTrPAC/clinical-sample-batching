#!/bin/R
# reformat example inputs 

library(data.table)
library(readxl)
library(writexl)

# reformat read manifests 
adu830 = fread("ADU830-10060.csv", sep=",")
stanford_adu830 = read_shipment("Stanford_ADU830-10106_062121.xlsx")
ped830 = fread("PED830-10062.csv", sep=",")
stanford_ped830 = read_shipment("Stanford_PED830-10107_062121.xlsx")

manifest = rbindlist(list(adu830, ped830), fill=T)
shipment = rbindlist(list(stanford_adu830, stanford_ped830), fill=T)

# remove unused column names 
manifest = manifest[,.(ManifestID, Study, codedSiteID, PID, VialLabel, barcode, visitCode, SampleTypeCode, randomGroupCode, sex_psca, calculatedAge)]
shipment = shipment[,.(Box, Position, Viallabel, `2D Barcode`, `ppt type`, `Sample Type`)]

## change PID, viallabel, barcode
# PID 
pids = unique(manifest[,PID])
new_pids = sample(1e7:2e7, size=length(pids), replace=FALSE)
names(new_pids) = pids
manifest[,PID := new_pids[as.character(PID)]]

# viallabel 
vl = c(unique(manifest[,VialLabel]), unique(shipment[,Viallabel]))
new_vl = as.numeric(sample(1e10:2e10, size=length(vl), replace=FALSE))
names(new_vl) = vl
manifest[,VialLabel := new_vl[as.character(VialLabel)]]
shipment[,Viallabel := new_vl[as.character(Viallabel)]]

# barcode 
bc = c(unique(manifest[,barcode]), unique(shipment[,`2D Barcode`]))
new_bc = as.numeric(sample(1e9:2e9, size=length(bc), replace=FALSE))
names(new_bc) = bc
manifest[,barcode := new_bc[as.character(barcode)]]
shipment[,`2D Barcode` := new_bc[as.character(`2D Barcode`)]]

# write out 
write.table(manifest, "examples/input/biospecimen-metadata.csv", sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE)
write_xlsx(shipment, "examples/input/shipment-manifest.xlsx")
