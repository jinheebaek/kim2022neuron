##
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

##
zbase = [-1.5, -1]
t_plot = [-1.5, 2]
figroot = "tfpower/around_freezing_onset/figure/grp_relz"

chsel = "pow_med"

##
tfr = pd.read_feather("tfpower/around_freezing_onset/data/morlet_wavg_rel_{}.feather".format(chsel))
freqs = tfr.freq.unique()
t1 = tfr.time.min()
t2 = tfr.time.max()


def norm(data):
    p_rel = data.pivot(index="freq", columns="time", values="power_rel")

    m_rel = p_rel.loc[:, zbase[0]:zbase[1]].mean(1).values.reshape(-1, 1)
    s_rel = np.std(p_rel.loc[:, zbase[0]:zbase[1]], axis=1, ddof=1).values.reshape(-1, 1)
    p_rel_z = (p_rel - m_rel) / s_rel

    p_rel_z = p_rel_z.stack().reset_index().rename(columns={0: "power_rel_z"})

    result = (
        p_rel_z
        # .merge(p_rel, on=["freq", "time"])
        # .merge(p_rel_z, on=["freq", "time"])
        .reset_index(drop=True)
    )
    return result


tfr_sbj = (
    tfr.groupby(["experiment", "region", "subject", "session"])
    .apply(norm).reset_index()
    .merge(tfr, on=["experiment", "region", "subject", "session", "freq", "time"])
    .loc[:, ["experiment", "region", "subject", "session", "freq", "time",
             "power_rel", "power_rel_z"]]
)

tfr_grp = (
    tfr_sbj
    .groupby(["experiment", "region", "session", "freq", "time"])
    .mean().reset_index()
)

del tfr
del tfr_sbj

##
fig, axes = plt.subplots(1, 2, figsize=(14, 4), constrained_layout=True)
side2axidx = {
    "l": 0, "r": 1
}
for (exp, ses), sbjdata in tfr_grp.groupby(["experiment", "session"]):
    print("processing exp: {}, ses: {}".format(exp, ses))

    for reg, regdata in sbjdata.groupby("region"):
        side = reg[0]
        sideidx = side2axidx[side]

        if reg.endswith("acc"):
            vmin = -8
            vmax = 11
        elif reg.endswith("bla"):
            vmin = -6
            vmax = 7

        p_rel_z = regdata.pivot(index="freq", columns="time", values="power_rel_z")
        ax = axes[sideidx]
        p_img = ax.imshow(
            p_rel_z, cmap="jet",
            vmin=vmin, vmax=vmax,
            origin="lower", aspect="auto",
            extent=np.array([t1, t2, freqs[0], freqs[-1]])
        )
        del p_rel_z
        del regdata

    for _, ax in np.ndenumerate(axes):
        ax.set_ylim([freqs.min(), 16])
        ax.set_yticks([4, 8, 12, 16])
        ax.set_yticklabels(ax.get_yticks(), fontsize=32, family="Arial")

        ax.set_xlim(t_plot)
        ax.set_xticks(np.arange(t_plot[0], t_plot[-1] + 0.5, 0.5))
        ax.set_xticklabels(ax.get_xticks(), fontsize=32, family="Arial")

        ax.set_ylabel("Frequency (Hz)", fontsize=36, family="Arial")
        ax.set_xlabel("Time from freezing onset (s)", fontsize=36, family="Arial")

        ax.axvline(0, ls="--", color="black", lw=3)
        ax.axhline(3, ls="--", color="white", lw=3)
        ax.axhline(5, ls="--", color="white", lw=3)
        ax.axhline(7, ls="--", color="white", lw=3)

    figdir = "{}/{}".format(chsel, exp)
    figdir = os.path.join(figroot, figdir)

    if not os.path.exists(figdir):
        os.makedirs(figdir)

    # figfile = "{}.png".format(ses)
    figfile = "{}.svg".format(ses)
    fig.savefig(
        os.path.join(figdir, figfile), transparent=True
    )

    for _, ax in np.ndenumerate(axes):
        ax.cla()

    del sbjdata

    fig_cbar, ax_cbar = plt.subplots(1, 1, figsize=(3, 10), dpi=300, facecolor="none")
    ax_cbar.axis("off")
    cbar = plt.colorbar(
        # p_img, ax=ax_cbar, ticks=np.arange(-3, 5, 1), aspect=30#, pad=0.5, anchor=(0.5, 0.5)
        p_img, ax=ax_cbar, ticks=[vmin, 0, vmax], aspect=30#, pad=0.5, anchor=(0.5, 0.5)
    )
    cbar.ax.set_yticklabels(
        cbar.ax.get_yticks(),
        fontsize=24, family="Arial")
    fig_cbar.tight_layout()
    fig_cbar.savefig(
        os.path.join(figdir, "cbar.svg"), transparent=True
    )
    plt.close(fig_cbar)

plt.close()
print("done")
