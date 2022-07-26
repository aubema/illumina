#!/usr/bin/env python3

import numpy as np
import pandas as pd
import yaml

with open("iss_params.in") as f:
    p = yaml.safe_load(f)

df = pd.read_csv(f"{p['wd']}/obs.csv")

mask = (df["distance"] >= p["distance_min"]) & (
    df["distance"] <= p["distance_max"]
)

out = dict(
    longitude=df["lons"][mask].to_list(),
    latitude=df["lats"][mask].to_list(),
)

with open("coordinates.dat", "w") as f:
    f.write(yaml.safe_dump(out))
