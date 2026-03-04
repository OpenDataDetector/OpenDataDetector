#!/usr/bin/env python3

import matplotlib.pyplot as plt
import mplhep
import argparse
from pathlib import Path
import json
import math
import ROOT

from root_helper import TH1


plt.rcParams["ytick.right"] = plt.rcParams["xtick.top"] = True
plt.rcParams["xtick.direction"] = "in"
plt.rcParams["ytick.direction"] = "in"
plt.rcParams["font.size"] = 12.0
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["legend.frameon"] = False
plt.rcParams["legend.columnspacing"] = 0.2
plt.rcParams["legend.handletextpad"] = 0.2
plt.rcParams["legend.labelspacing"] = 0.2
plt.rcParams["legend.borderpad"] = 0
plt.rcParams["legend.handlelength"] = 1.0


parser = argparse.ArgumentParser(description="Make material composition plots")
parser.add_argument("input", type=Path, help="Input root file with histograms")
parser.add_argument("config", type=Path, help="Input config file with colours")
parser.add_argument("output", type=Path, help="Output directory for plots")
args = parser.parse_args()

args.output.mkdir(parents=True, exist_ok=True)

rf = ROOT.TFile.Open(args.input.absolute().as_posix())

tracker = {
    "beampipe": "Beam pipe",
    "pixel": "Pixel",
    "sstrips": "Short Strips",
    "lstrips": "Long Strips",
    "solenoid": "Solenoid",
}

calo = {
    "ecalbarrel": "ECal barrel",
    "ecalendcap": "ECal endcap",
    "hcalbarrel": "HCal barrel",
    "hcalendcap": "HCal endcap",
}

with args.config.open() as f:
    config = json.load(f)
configcolours = {det: config[det][0]["colour"] for det in config.keys()}

for group, names in [
    ("tracker_", {**tracker}),
    ("calo_", {**calo}),
    ("", {**tracker, **calo}),
]:
    for y in ("x0", "l0"):
        for x in ("phi", "eta"):
            x_lim = {"phi": (-math.pi, math.pi), "eta": (-4, 4)}[x]

            hists = []
            bins = None
            labels = []
            colours = []
            for name, label in names.items():
                key = f"{name}_{y}_vs_{x}_all"
                if rf.Get(key) is not None:
                    th1 = TH1(rf.Get(key), xrange=x_lim)
                    edges = list(th1.x_lo) + [th1.x_hi[-1]]

                    hists.append(th1.y)
                    bins = edges
                    labels.append(label)
                    colours.append(configcolours[name])
                else:
                    print(f"Key {key} not found in {args.input}")

            fig, ax = plt.subplots()
            mplhep.histplot(
                hists,
                bins=bins,
                ax=ax,
                stack=True,
                histtype="fill",
                label=labels,
                color=colours,
            )
            ymin, ymax = ax.get_ylim()
            ax.set_ylim(top=1.2 * ymax)
            ax.legend(ncol=3)

            ylab = {"l0": r"\lambda_0", "x0": "X_0"}[y]
            ax.set_ylabel(rf"${ylab}$")
            ax.set_xlim(bins[0], bins[-1])
            fig.savefig(args.output / f"{group}{y}_vs_{x}.pdf")
