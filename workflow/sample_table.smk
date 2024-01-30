pepfile: "sample_table_config.yaml"


SAMPLES = pep.sample_table["sample_name"]
PAIRED = pep.sample_table.columns.str.contains("R2").any()

if PAIRED:
    FRACTIONS = ["R1", "R2"]
else:
    FRACTIONS = ["se"]


def get_quality_controlled_reads(wildcards):
    headers = ["Reads_QC_" + f for f in FRACTIONS]
    
    return pep.sample_table.loc[wildcards.sample, headers].tolist()


# def get_quality_controlled_reads(wildcards):
#     return (
#         expand(
#             "QC/reads/{sample}_{fraction}.fastq.gz",
#             fraction=FRACTIONS,
#             sample=wildcards.sample,
#         ),
#     )
