ref_fa="test/data/"
query_fa="test/data/"
output_directory="$(dirname $ref_fa)/results"

echo -e "======\n Testing NF execution \n======" \
&& rm -rf $output_directory \
&& nextflow run nf-align_proteins.nf \
	--ref_fa $ref_fa \
  --query_fa $query_fa \
	--output_dir $output_directory \
	-resume \
	-with-report $output_directory/`date +%Y%m%d_%H%M%S`_report.html \
	-with-dag $output_directory/`date +%Y%m%d_%H%M%S`.DAG.html \
	-with-timeline $output_directory/`date +%Y%m%d_%H%M%S`_timeline.html \
&& echo -e "======\n Basic pipeline TEST SUCCESSFUL \n======"
	#-stub-run \
