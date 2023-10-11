#!/usr/bin/env python3

import uproot
import matplotlib.pyplot as plt
import numpy as np
import os
import hist
import mplhep
import argparse
from pathlib import Path

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
parser.add_argument("input", help="Input root file with histograms")
parser.add_argument("output", type=Path, help="Output directory for plots")
args = parser.parse_args()

args.output.mkdir(parents=True, exist_ok=True)

rf = uproot.open(args.input)

names = {
    "beampipe": "Beam pipe",
    "sstrips": "Short Strips",
    "lstrips": "Long Strips",
    "pixel": "Pixel",
    "solenoid": "Solenoid",
    "ecalbarrel": "EMCal barrel",
    "ecalrightendcap": "EMCal right endcap",
    "ecalleftendcap": "EMCal left endcap"
}

for q in (
    "x0", 
    "l0",
):
    for qq in ("phi", "eta"):
        h = None
        hists = []
        labels = []
        detector = None
        for k in rf:
            name, _ = k.split(";", 1)
            if not name.endswith("all"): continue
            if not q in name: continue
            if not qq in name: continue
            if name.startswith("detector"):
                detector = values.copy()
                continue

            o = rf[k].to_hist()
            o.axes[0].label = qq
            hists.append(o)
            l, _ = k.split("_", 1)
            l = names[l]
            labels.append(l)
            
        fig, ax = plt.subplots()
        mplhep.histplot(hists, ax=ax, stack=True, histtype="fill", label=labels)
        ymin, ymax = ax.get_ylim()
        ax.set_ylim(top=1.2*ymax)
        ax.legend(ncol=3)

        ylab = {"l0": "\lambda_0", "x0": "X_0"}[q]
        ax.set_ylabel(rf"${ylab}$")
        ax.set_xlim(hists[0].axes[0].edges[0], hists[0].axes[0].edges[-1])
        
        fig.savefig(args.output / f"{q}_vs_{qq}.pdf")

