##
import numpy as np
import pandas as pd
import math
import multiprocessing as mp

##
n_cycles = 5
fs = 1000 # approx.

##
# (t)  -->  (floor((t - window_size + 1) / step_size), window_size))
def rolling_window(a, window_size, step_size):
    shape = (int((a.shape[0] - window_size + 1) / step_size), window_size)
    strides = (a.strides[0] * step_size, a.strides[0])
    return np.lib.stride_tricks.as_strided(
        a, shape=shape, strides=strides, writeable=False)


def plv_over_time_each_freq(freqdata):
    freq = freqdata.index.get_level_values("freq").unique()
    assert np.size(freq) == 1
    freq = freq[0]

    window_size = math.ceil(1 / freq * n_cycles * fs)
    if window_size % 2 == 0:
        window_size += 1
    step_size = 1

    windowed = rolling_window(freqdata.values, window_size, step_size)

    plv = np.apply_along_axis(
        lambda x: np.abs(np.mean(np.exp(1j * np.angle(x)))),
        1, windowed
    )
    # t = rolling_window(freqdata.index.get_level_values("time").values,
    #                    window_size, step_size).mean(1)
    t = rolling_window(freqdata.index.get_level_values("time").values,
                       window_size, step_size)[:, math.ceil(window_size / 2 - 1)]
    return pd.Series(
        plv, name="plv",
        index=pd.Index(t, name="time")
    )


def calc_plv(dfseg):
    sbj = dfseg.index.get_level_values("subject").unique()
    assert np.size(sbj) == 1
    sbj = sbj[0]

    iepc = dfseg.index.get_level_values("iepoch").unique()
    assert np.size(iepc) == 1
    iepc = iepc[0]

    onset = dfseg.index.get_level_values("onset").unique()
    assert np.size(onset) == 1
    onset = onset[0]

    freezdur = dfseg.index.get_level_values("freezdur").unique()
    assert np.size(freezdur) == 1
    freezdur = freezdur[0]

    session = dfseg.index.get_level_values("session").unique()
    assert np.size(session) == 1
    session = session[0]

    plvdata = (
        dfseg
        .groupby("freq")
        .apply(plv_over_time_each_freq)
        .to_frame()
    )
    plvdata["subject"] = sbj
    plvdata["iepoch"] = iepc
    plvdata["onset"] = onset
    plvdata["session"] = session
    plvdata["freezdur"] = freezdur

    return plvdata


##
# plv = []
cwt_data_nodes = pd.read_csv("phase/around_freezing_onset/data/cwt_nodes.csv")
for _, node in cwt_data_nodes.iterrows():
    session = node.session
    sbj = "{:02.0f}".format(node.subject)

    if session != "observational_fear":
        continue
    # if sbj not in ("07", "22", "27"):
    #     continue

    print("processing session: {}, subject: {}".format(session, sbj))
    tfr_morlet = pd.read_hdf(node.data_file, node.data_node)

    data_reg1 = (
        tfr_morlet
        .query("region == '{}'".format("racc"))
        .set_index(["session", "subject", "iepoch", "onset", "freezdur", "time", "freq"])
    )
    data_reg2 = (
        tfr_morlet
        .query("region == '{}'".format("rbla"))
        .set_index(["session", "subject", "iepoch", "onset", "freezdur", "time", "freq"])
    )

    csd = data_reg1.sig * np.conjugate(data_reg2.sig)

    with mp.Pool(mp.cpu_count()) as p:
        plvdata = p.map(calc_plv,
                        [grp for (sbj, iepch), grp in csd.groupby(
                            ["subject", "iepoch"])])

    plvdata = (
        pd.concat(plvdata)
        # .groupby(["batch", "subject", "freq", "time"])
        # .mean()                 # average epochs
        # .reset_index()
    )

    plvdata["reg1"] = "racc"
    plvdata["reg2"] = "rbla"
    plvdata = plvdata.rename(columns={"sig": "plv"})
    # plv.append(plvdata)

    ## to hdf
    session = node.session
    subject = "{:02.0f}".format(node.subject)
    plvdata.reset_index().to_hdf(
        "phase/around_freezing_onset/data/plv.h5", "{}/{}".format(
            session, subject)
    )

# plv = pd.concat(plv).reset_index()
# plv.to_feather("phase/around_freezing_onset/data/plv.feather")

print("done")
