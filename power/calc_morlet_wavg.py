##
import os
import neuro
import mne

import numpy as np
import pandas as pd

##
t_window = [-2.1, 2.1]

freqs = np.arange(2, 20, 0.2)
ncycles = 5

##
data_dir = "tfpower/around_freezing_onset/data"
data_file = "morlet_wavg_rel_{}.feather"
chsel = "pow_med"
data_file = data_file.format(chsel)

##
lfp_data_nodes = pd.read_csv("lfp/data/data_nodes.csv")
lfp_ch_candidates = pd.read_csv("lfp/data/channel_candidates.csv")
excludes = pd.read_csv("freezing/data/freezing_w_noise_{}.csv".format(chsel))

freezing = (
    pd.read_feather("freezing/data/freezing_all.feather")
    .assign(
        behab_session=lambda x: np.where(
            x.session == "24h_memory", "24h", np.where(
            x.onset < 180, "hab", "of"
        )),
        t_from_prev_offset=lambda x: x.onset - x.offset.shift(1),
        t_from_prev_onset=lambda x: x.onset - x.onset.shift(1)
    )
    ## exclude theta with interval <1.5 sec (baseline period)
    .query("(t_from_prev_offset < 0) | (t_from_prev_offset >= 1.5) | (t_from_prev_offset != t_from_prev_offset)")
)

for _, exepc in excludes.iterrows():
    freezing = (
        freezing
        ## exclude noisy segs
        .query("~((experiment==@exepc.experiment)&(subject==@exepc.subject)&(session==@exepc.session)&(onset==@exepc.onset))")
    )

lfp_data_nodes = (
    lfp_data_nodes
    # .query("session == 'observational_fear'")
    .merge(lfp_ch_candidates, how="left", on=["experiment", "session", "subject", "region"])
)

##
tfr_morlet_wavg = []
for (exp, sbj, behabses), frzs in freezing.groupby(["experiment", "subject", "behab_session"]):
    print("processing exp: {}, sbj: {}, ses: {}".format(exp, sbj, behabses))
    ses = frzs.session.unique()
    assert ses.shape[0] == 1
    ses = ses[0]

    lfp_nodes = (
        lfp_data_nodes
        .query("experiment == @exp")
        .query("subject == @sbj")
        .query("session == @ses")
    )

    lfps = []
    for _, node in lfp_nodes.iterrows():
        lfp_data = neuro.read_stream(node.data_file, node.data_node)
        lfps.append({
            "region": node.region,
            "lfp_data": lfp_data._data.loc[:, node["ch_by_" + chsel]]
        })
        fs = lfp_data.fs
        del lfp_data
    lfps = pd.DataFrame(lfps)

    nchannels = lfp_nodes.shape[0]
    ntimes = int((t_window[1] - t_window[0]) * fs)
    nepoch = frzs.shape[0]

    epoch_data_onset = np.empty((nepoch, nchannels, ntimes))

    iepoch = 0
    for _, frz_data in frzs.iterrows():
        t_onset = frz_data.onset

        ## lfp
        ireg = 0
        for _, regdata in lfps.iterrows():
            lfp_seg_tmp = regdata.lfp_data
            lfp_seg = lfp_seg_tmp[(t_onset + t_window[0]):]
            lfp_seg = lfp_seg.iloc[:ntimes]
            del lfp_seg_tmp

            t = np.arange(ntimes) * 1/fs
            idxt0 = int(abs(t_window[0]) * fs)
            t = t - t[idxt0]
            epoch_data_onset[iepoch, ireg, :] = lfp_seg

            ireg += 1
        iepoch += 1

    tfr_onset_epoch = mne.time_frequency.tfr_array_morlet(
        epoch_data_onset, float(fs), freqs, n_jobs=6, output="power", n_cycles=ncycles
    )

    for iepch in range(tfr_onset_epoch.shape[0]):
        for ich in range(tfr_onset_epoch.shape[1]):
            tfr_onset_epoch[iepch, ich, :, :] /= tfr_onset_epoch[iepch, ich, :, :].sum(0)
    tfr_onset_epoch *= 100

    ## weighted average
    tfr_onset_wavg = np.average(
        tfr_onset_epoch, axis=0,
        weights=frzs.duration
    )

    for ich in range(tfr_onset_wavg.shape[0]):
        tfr_onset_wavg_ch = tfr_onset_wavg[ich]
        tfr_onset_wavg_ch = pd.DataFrame(
            tfr_onset_wavg_ch,
            index=pd.Series(freqs, name="freq"),
            columns=pd.Series(t, name="time")
        )
        tfr_onset_wavg_ch = (
            tfr_onset_wavg_ch
            .stack().reset_index()
            .rename(columns={0: "power_rel"})
        )
        tfr_onset_wavg_ch["experiment"] = exp
        tfr_onset_wavg_ch["region"] = lfps.region.iloc[ich]
        tfr_onset_wavg_ch["subject"] = sbj
        tfr_onset_wavg_ch["session"] = behabses
        tfr_onset_wavg_ch["t0"] = "onset"
        tfr_morlet_wavg.append(tfr_onset_wavg_ch)
        del tfr_onset_wavg_ch

    del tfr_onset_epoch
    del lfps

tfr_morlet_wavg = pd.concat(tfr_morlet_wavg).reset_index(drop=True)

if not os.path.exists(data_dir):
    os.makedirs(data_dir)

tfr_morlet_wavg.reset_index(drop=True).to_feather(
    os.path.join(data_dir, data_file.format(chsel))
    # , compression="uncompressed"
)

print("done")
