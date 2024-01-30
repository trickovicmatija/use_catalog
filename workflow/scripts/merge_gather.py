
import os, sys
import logging, traceback

logging.basicConfig(
    filename=snakemake.log[0],
    level=logging.DEBUG,
    format="%(asctime)s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

logging.captureWarnings(True)


def handle_exception(exc_type, exc_value, exc_traceback):
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return

    logging.error(
        "".join(
            [
                "Uncaught exception: ",
                *traceback.format_exception(exc_type, exc_value, exc_traceback),
            ]
        )
    )


# Install exception handler
sys.excepthook = handle_exception

### Start of script


import pandas as pd
from pathlib import Path
import os

import sys

data_file=Path(snakemake.output[0])


input_files = list(snakemake.input)


combined = dict()


for sample_path in input_files:

    sample_file = os.path.basename(sample_path)
    sample = os.path.splitext(sample_file)[0]

    df = pd.read_csv(sample_path)

    #df = df.loc[df.potential_false_negative == False] # not available in my version
    df = df.loc[df.intersect_bp >= 25000]
    #df = df.loc[df['f_match'] > 0.01] # I experiment with there, but at the end keep it commeted out

    combined[sample]=df.set_index("name").median_abund


D = pd.concat(combined,axis=1).fillna(0).astype(int)

D.index= D.index.astype(str).str.rjust(10,"0") # add leading zeros to subspecies names

D.T.to_hdf(data_file,"median_abund")