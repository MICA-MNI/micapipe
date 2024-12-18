#!/usr/bin/env python3
from pathlib import Path

from snakebids import bidsapp, plugins


app = bidsapp.app(
    [
        plugins.SnakemakeBidsApp(Path(__file__).resolve().parent),
        plugins.BidsValidator(),
        plugins.Version(distribution="micapipe_snake"),
        plugins.CliConfig("parse_args"),
        plugins.ComponentEdit("pybids_inputs"),
    ]
)


if __name__ == "__main__":
    app.run()