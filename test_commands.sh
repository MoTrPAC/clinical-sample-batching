Rscript randomization.R \
    --shipment-manifest-excel \
    	/Users/nicolegay/Documents/motrpac/CLINICAL/batching_officer/Stanford_ADU830-10060_120720.xlsx \
    	/Users/nicolegay/Documents/motrpac/CLINICAL/batching_officer/Stanford_PED830-10062_120720.xlsx \
    --api-metadata-csv \
    	/Users/nicolegay/Documents/motrpac/CLINICAL/batching_officer/ADU830-10060.csv \
    	/Users/nicolegay/Documents/motrpac/CLINICAL/batching_officer/PED830-10062.csv \
    --max-n-per-batch 94 \
    --outdir ~/Desktop/stanford_batches 

Rscript randomization.R \
    -ship \
        /Users/nicolegay/Documents/motrpac/CLINICAL/batching_officer/Stanford_ADU830-10060_120720.xlsx \
        /Users/nicolegay/Documents/motrpac/CLINICAL/batching_officer/Stanford_PED830-10062_120720.xlsx \
    -api \
        /Users/nicolegay/Documents/motrpac/CLINICAL/batching_officer/ADU830-10060.csv \
        /Users/nicolegay/Documents/motrpac/CLINICAL/batching_officer/PED830-10062.csv \
    -max 94 \
    -o ../batches > ~/Desktop/stanford_batches/out.log 2>&1

Rscript randomization.R \
    -ship ~/Desktop/broad_batches/ShipmentContents_BroadCarr_012521.xlsx \
    -api ~/Desktop/broad_batches/ADU822-10074.csv \
    -max 15 \
    -s \
    -o ~/Desktop/broad_batches > ~/Desktop/broad_batches/out.log 2>&1