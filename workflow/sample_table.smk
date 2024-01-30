import pandas as pd

ADDITIONAL_SAMPLEFILE_HEADERS = []


def validate_sample_table(sampleTable):

    Expected_Headers = ADDITIONAL_SAMPLEFILE_HEADERS
    for h in Expected_Headers:
        if not (h in sampleTable.columns):
            Exception(f"expect '{h}' to be found in samples.tsv")

    if not sampleTable.index.is_unique:
        duplicated_samples = ", ".join(D.index.duplicated())
        Exception(
            f"Expect Samples to be unique. Found {duplicated_samples} more than once"
        )

    if not (sampleTable.columns.str.startswith("Reads_QC_").sum() >= 1):

        raise IOError(
            "QC reads need to be in the sample table. "
            "I din't find any colums with Reads_QC_<fraction>'  "
        )


def load_sample_table(sample_table="samples.tsv"):

    sampleTable = pd.read_csv(sample_table, index_col=0, sep="\t")
    validate_sample_table(sampleTable)
    return sampleTable


sampleTable = load_sample_table(config.get("sample_table", "samples.tsv"))

SAMPLES = sampleTable.index.values
if sampleTable.columns.str.contains("R2").any():
    MULTIFILE_FRACTIONS = ["R1", "R2"]
else:
    MULTIFILE_FRACTIONS = ["se"]


class FileNotInSampleTableException(Exception):
    """
    Exception with sampleTable
    """

    def __init__(self, message):
        super(FileNotInSampleTableException, self).__init__(message)


def get_files_from_sampleTable(sample, Headers):
    """
    Function that gets some filenames form the sampleTable for a given sample and Headers.
    It checks various possibilities for errors and throws either a
    FileNotInSampleTableException or a IOError, when something went really wrong.
    """

    if not (sample in sampleTable.index):
        raise FileNotInSampleTableException(
            f"Sample name {sample} is not in sampleTable"
        )

    Error_details = f"\nsample: {sample}\nFiles: {Headers}"

    if type(Headers) == str:
        Headers = [Headers]

    NheadersFound = sampleTable.columns.isin(Headers).sum()

    if NheadersFound == 0:
        raise FileNotInSampleTableException(
            f"None of the Files ar in sampleTable, they should be added to the sampleTable later in the workflow"
            + Error_details
        )
    elif NheadersFound < len(Headers):
        raise IOError(
            f"Not all of the Headers are in sampleTable, found only {NheadersFound}, something went wrong."
            + Error_details
        )

    files = sampleTable.loc[sample, Headers]

    if files.isnull().all():
        raise FileNotInSampleTableException(
            "The following files were not available for this sample in the SampleTable"
            + Error_details
        )

    elif files.isnull().any():
        raise IOError(
            f"Not all of the files are in sampleTable, something went wrong."
            + Error_details
        )

    return list(files)


def get_raw_reads(wildcards):
    """Gets quality controlled reads. Two files if paired end and only one if single end."""

    Headers = ["Reads_raw_" + f for f in MULTIFILE_FRACTIONS]
    sample_dir_path = "/".join(config["sample_table"].split("/")[:-1]) + "/"
    return [
        sample_dir_path + s
        for s in get_files_from_sampleTable(wildcards.sample, Headers)
    ]



def get_quality_controlled_reads(wildcards):
    """Gets quality controlled reads. Two files if paired end and only one if single end."""

    QC_Headers = ["Reads_QC_" + f for f in MULTIFILE_FRACTIONS]
    sample_dir_path = "/".join(config["sample_table"].split("/")[:-1]) + "/"
    return [
        sample_dir_path + s
        for s in get_files_from_sampleTable(wildcards.sample, QC_Headers)
    ]


