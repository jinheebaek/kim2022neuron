##
import os
import numpy as np
import pandas as pd

import scipy.stats
import matplotlib.pyplot as plt

##
t1 = 1.5
t2 = 2

base = (-1.5, -1)
##


def norm(data):
    p = data.pivot(index="freq", columns="time", values="plv")

    m = p.loc[:, base[0]:base[1]].mean(1).values.reshape(-1, 1)
    s = np.std(p.loc[:, base[0]:base[1]], axis=1, ddof=1).values.reshape(-1, 1)
    p_z = (p - m) / s

    p = p.stack().reset_index().rename(columns={0: "plv"})
    p_z = p_z.stack().reset_index().rename(columns={0: "plv_z"})
    return pd.merge(p_z, p, on=["freq", "time"]).reset_index(drop=True)


##
# delete prev. data before running
data_dst = "phase/around_freezing_onset/data/plv4stat.h5"
if os.path.exists(data_dst):
    os.unlink(data_dst)

#
with pd.HDFStore("phase/around_freezing_onset/data/plv.h5") as store:
    for k in store:
        meta = k.split("/")
        ses = meta[1]
        sbj = meta[2]

        if ses != "observational_fear":
            continue

        print("processing ses: {}, sbj: {}".format(ses, sbj))
        sbjdata = (
            store.get(k)
            .merge(
                (pd.read_csv("subject/group.csv")
                 .assign(subject=lambda x: x.subject.apply(lambda x: "{:02.0f}".format(x)))),
                how="left", on="subject")
        )
        grp = sbjdata.group.iloc[0]

        plv4stat = (
            sbjdata
            .groupby(["group", "session", "subject", "freq", "time"])
            .apply(lambda x: np.average(x.plv, weights=x.freezdur)).reset_index().rename(columns={0: "plv"})
            .pipe(norm)
            .assign(
                group=grp, session=ses, subject=sbj
            )
            .loc[:, ["group", "session", "subject", "freq", "time", "plv", "plv_z"]]
        )

        with pd.HDFStore("phase/around_freezing_onset/data/plv4stat.h5") as store_dst:
            store_dst.append("data", plv4stat,
                             min_itemsize={
                                 "group": len("post-theta"),
                                 "session": len("observational_fear"),
                                 "subject": 5
                             },
                             format="table", data_columns=True)

        del sbjdata
        del plv4stat

##
plv4plot = (
    pd.read_hdf("phase/around_freezing_onset/data/plv4stat.h5", "data/table")
    .reset_index()
    .groupby(["group", "session", "freq", "time"])
    .mean()
    .reset_index()
)
plv4plot.to_feather("phase/around_freezing_onset/data/plv4plot.feather")

print("done")
