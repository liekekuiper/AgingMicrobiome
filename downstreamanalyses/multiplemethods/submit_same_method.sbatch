#!/bin/bash

curl -O http://ftp.microbio.me/greengenes_release/2022.10/2022.10.taxonomy.asv.tsv.qza
curl -O http://ftp.microbio.me/greengenes_release/2022.10/2022.10.phylogeny.asv.nwk.qza

version=2022.10
taxonomy=${version}.taxonomy.asv.tsv.qza
tree=${version}.phylogeny.asv.nwk.qza

depth_16S=10000
depth_WGS=50000
metadata=MrOS_meta.txt
table_16s=./16S/renamed_feature.table.gg2-2022.10.qza
table_wgs=./MrOS_Shotgun_Sequencing/renamed_feature.table.gg2-2022.10.qza
label=MrOS
sbatch --export depth_16S=${depth_16S},depth_WGS=${depth_WGS},tree=${tree},taxonomy=${taxonomy},metadata=${metadata},table_16s=${table_16s},table_wgs=${table_wgs},label=${label} submit.sbatch

