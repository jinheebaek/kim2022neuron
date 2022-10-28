##
import neuro
import scipy
import numpy as np
import pandas as pd

import mne

##
light = (
    pd.read_feather("light/data/473.feather")
    .query("onset >= 180")
    .query("offset <= 600")
)
group = pd.read_csv("subject/group.csv")
excludes = (
    pd.read_csv("freezing/data/freezing_w_noise.csv")
    .assign(subject=lambda x: x.subject.apply(lambda x: "{:02.0f}".format(x)))
)
freezing_of = (
    pd.read_csv("freezing/data/freezing.csv")
    .assign(
        data_session="observational_fear",
        session=lambda x: np.where(x.onset < 180, "habituation", "observational_fear"),
        t_from_prev_offset=lambda x: x.onset - x.offset.shift(1),
        t_from_prev_onset=lambda x: x.onset - x.onset.shift(1)
    )
    .query("(t_from_prev_offset < 0) | (t_from_prev_offset >= 1.5) | (t_from_prev_offset != t_from_prev_offset)")
)
freezing = (
    freezing_of
    .assign(subject=lambda x: x.subject.apply(lambda x: "{:02.0f}".format(x)))
    .reset_index(drop=True)
)
for _, exepc in excludes.iterrows():
    freezing = (
        freezing
        ## noise
        .query("~((subject == @exepc.subject)&(session == @exepc.session)&(onset == @exepc.onset))")
    )

lfp_data_nodes = (
    pd.read_csv("lfp/data/data_nodes.csv")
    .merge(group, how="left", on="subject")
    .assign(subject=lambda x: x.subject.apply(lambda x: "{:02.0f}".format(x)))
)


##
def stim_mask(eeg_series, t_stims, t_pre=0.01, t_post=0.05):
    "remove ChR2 artifact and interpolate date"
    x_new = eeg_series.index.values

    x = x_new
    for t in t_stims:
        x = x[((x < (t - t_pre)) | ((t + t_post) < x))]

    f = scipy.interpolate.interp1d(x, eeg_series.loc[x].values, kind="cubic")
    masked = pd.Series(f(x_new), index=x_new)

    return masked


##
t1 = 3
t2 = 3.5
freqs = np.arange(2, 20, 0.2)

cwt_nodes = []
for (sbj, data_session, session), frz_sbjses in freezing.groupby(["subject", "data_session", "session"]):
    if session != "observational_fear":
        continue

    print("processing session: {}, subject: {}".format(session, sbj))

    lfp_nodes = (
        lfp_data_nodes
        .query("subject == @sbj")
        .query("session == @data_session")
    )
    light_sbj = light.query("subject == @sbj")

    lfp_sbj = []
    assert(lfp_nodes.shape[0] == 1)
    data_node = lfp_nodes.iloc[0]
    lfp_data = neuro.read_stream(data_node.data_file, data_node.data_node)

    broad_raw = lfp_data.bandpass(1, 40)._data
    for reg, regdata in broad_raw.iteritems():
        if session == "observational_fear":
            broad_masked = stim_mask(broad_raw.loc[:, reg], light_sbj.onset.values)
        else:
            broad_masked = broad_raw.loc[:, reg]

        lfp_sbj.append({
            "region": reg,
            "lfp_data": broad_masked
        })
    lfp_sbj = pd.DataFrame(lfp_sbj)

    fs = lfp_data.fs

    nchannels = lfp_sbj.shape[0]
    ntimes = int((t1 + t2) * fs)
    nepoch = frz_sbjses.shape[0]

    epoch_data_onset = np.empty((nepoch, nchannels, ntimes))
    t = np.arange((t1 + t2) * fs) / fs - t1

    iepoch = 0
    for _, frz_data in frz_sbjses.iterrows():
        t_onset = frz_data.onset

        ## lfp
        ireg = 0
        for _, regdata in lfp_sbj.iterrows():
            lfp_data = regdata.lfp_data
            lfp_seg = lfp_data[(t_onset - t1):]
            lfp_seg = lfp_seg.iloc[:ntimes]

            # t = lfp_seg.index - t_onset
            t = np.arange(ntimes) * 1/fs
            idxt0 = int(abs(t1) * fs)
            t = t - t[idxt0]
            epoch_data_onset[iepoch, ireg, :] = lfp_seg

            del regdata
            del lfp_data
            del lfp_seg

            ireg += 1
        iepoch += 1

    tfr_sbj_onset = mne.time_frequency.tfr_array_morlet(
        epoch_data_onset, float(fs), freqs, n_jobs=5, output="complex", n_cycles=5
    )

    tfr_morlet_sbjses = []
    for iepoch in range(tfr_sbj_onset.shape[0]):
        for ich in range(tfr_sbj_onset.shape[1]):
            tfr_sbj_ch_onset = tfr_sbj_onset[iepoch, ich, :, :]
            tfr_sbj_ch_onset = pd.DataFrame(
                tfr_sbj_ch_onset,
                index=pd.Series(freqs, name="freq"),
                columns=pd.Series(t, name="time")
            )
            tfr_sbj_ch_onset = (
                tfr_sbj_ch_onset
                .stack().reset_index()
                .rename(columns={0: "sig"})
            )
            tfr_sbj_ch_onset["session"] = session
            tfr_sbj_ch_onset["region"] = lfp_sbj.region.iloc[ich]
            tfr_sbj_ch_onset["subject"] = sbj
            tfr_sbj_ch_onset["t0"] = "onset"
            tfr_sbj_ch_onset["iepoch"] = iepoch
            tfr_sbj_ch_onset["onset"] = frz_sbjses.onset.iloc[iepoch]
            tfr_sbj_ch_onset["freezdur"] = frz_sbjses.duration.iloc[iepoch]
            # tfr_morlet.append(tfr_sbj_ch_onset)
            tfr_morlet_sbjses.append(tfr_sbj_ch_onset)

            del tfr_sbj_ch_onset

    tfr_morlet_sbjses = pd.concat(tfr_morlet_sbjses).reset_index(drop=True)
    tfr_morlet_sbjses.to_hdf("phase/around_freezing_onset/data/cwt.h5", "{}/{}".format(sbj, session))
    del tfr_morlet_sbjses

    cwt_nodes.append({
        "data_file": "phase/around_freezing_onset/data/cwt.h5",
        "data_node": "{}/{}".format(sbj, session),
        "subject": sbj,
        "session": session
    })

    del lfp_sbj
    del tfr_sbj_onset
    del epoch_data_onset

# tfr_morlet = pd.concat(tfr_morlet).reset_index(drop=True)
# tfr_morlet.to_hdf("phase/around_freezing_onset/data/cwt.h5", "data")

cwt_nodes = pd.DataFrame(cwt_nodes)
cwt_nodes.to_csv("phase/around_freezing_onset/data/cwt_nodes.csv", index=False)

print("done")
