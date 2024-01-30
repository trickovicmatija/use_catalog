# use_catalog
Use specific sourmash index to quantify bacterial taxa from metagenomic samples.


The pipeline uses [PEP specification](https://pep.databio.org/en/latest/) to specify the file paths.

To get started.
Copy the `config/sample_table_config.yaml`  and file to your working directory.
Create a sample_table.csv for your samples. similar to the [`config/sample_table.csv`](config/sample_table.csv) file. 

The first column needs to be `sample_name`. I'am looking for the headers `Reads_QC_R1`, `Reads_QC_R2` or `Reads_QC_se` 

If you have the fastq files under a different header you can rename the columns using the PEP specification. 



