##
import numpy as np
import pandas as pd

import scipy.stats
import matplotlib.pyplot as plt

##
t1 = 1.5
t2 = 2

base = (-1.5, -1)

vmin = -6
vmax = 4

##
plv4plot = pd.read_feather("phase/around_freezing_onset/data/plv4plot.feather")
freqs = plv4plot.freq.unique()

##
fig, ax = plt.subplots(1, 1, figsize=(8, 6), constrained_layout=True)
for grp, grpdata in plv4plot.query("session == 'observational_fear'").groupby("group"):
    # p = grpdata.pivot(index="freq", columns="time", values="plv")
    p = grpdata.pivot(index="freq", columns="time", values="plv_z")
    p = p.loc[:, -t1:t2]

    p_img = ax.imshow(
        p.values, aspect="auto",
        # cmap="jet",
        origin="lower",
        vmin=vmin, vmax=vmax,
        extent=np.array([-t1, t2, freqs[0], freqs[-1]])
    )

    ax.set_ylim([2, 16])
    ax.set_yticks([4, 8, 12, 16])
    ax.set_yticklabels(ax.get_yticks(), fontsize=38, family="Arial")

    # ax.set_xticks(np.arange(-t1, t2 + 1))
    ax.set_xticks([-1, 0, 1, 2])
    ax.set_xticklabels(ax.get_xticks(), fontsize=38, family="Arial")

    ax.set_ylabel(
        "Frequency (Hz)", fontsize=46, family="Arial")

    ax.set_xlabel(
        "Time from freezing onset (s)", fontsize=46, family="Arial"
    )
    # ax.set_title("{}".format(ses), fontsize=40)

    ax.axvline(0, ls="--", color="black", lw=3)
    ax.axhline(3, ls="--", color="white", lw=3)
    ax.axhline(5, ls="--", color="white", lw=3)
    ax.axhline(7, ls="--", color="white", lw=3)

    fig.savefig(
        "phase/around_freezing_onset/figure/plv_over_time_{}.svg".format(grp),
        transparent=True)

    ax.cla()
    # cbar.remove()
    del p

    fig_cbar, ax_cbar = plt.subplots(1, 1, figsize=(2, 3), dpi=300, facecolor="none", constrained_layout=True)
    ax_cbar.axis("off")
    cbar = plt.colorbar(
        p_img, ax=ax_cbar, ticks=[vmin, 0, vmax], aspect=15#, pad=0.5, anchor=(0.5, 0.5)
    )
    cbar.ax.set_yticklabels(
        cbar.ax.get_yticks(),
        fontsize=20, family="Arial")
    cbar.ax.set_ylabel("Phase locking value\n(z-score)", fontsize=18, labelpad=30, family="Arial", rotation=270)

    # fig_cbar.tight_layout()
    fig_cbar.savefig(
        "phase/around_freezing_onset/figure/cbar_{}.svg".format(grp),
        transparent=True
    )
    plt.close(fig_cbar)

plt.close(fig)

print("done")
