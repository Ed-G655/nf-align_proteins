##!/usr/bin/env bash
cd ../

mirna_gff="/bodega/projects/miRNome_project/references/gff/hsa.gff3"
utrbed="/bodega/projects/miRNome_project/references/beds/3utr.bed"
vcf="/bodega/projects/miRNome_project/references/vcfs/76g_PASS_0_05.vcf.gz"
fasta="/bodega/projects/miRNome_project/references/fastas/"
output_directory="/bodega/projects/miRNome_project/analisis/100GMX/results"

echo -e "======\n Lauching NF execution \n======" \
&& nextflow run nf-compare-miRNome.nf \
        --mirna_gff $mirna_gff \
        --utrbed $utrbed \
        --vcf $vcf \
        --fasta_dir $fasta \
        --output_dir $output_directory \
        -resume \
        -with-report $output_directory/`date +%Y%m%d_%H%M%S`_report.html \
        -with-dag $output_directory/`date +%Y%m%d_%H%M%S`.DAG.html \
        -with-timeline $output_directory/`date +%Y%m%d_%H%M%S`_timeline.html \
&& echo -e "======\n  pipeline execution END \n======"
