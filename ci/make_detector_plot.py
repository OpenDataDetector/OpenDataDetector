#!/usr/bin/env python3

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
import argparse
import json
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
parser.add_argument("config", type=Path, help="Input config file with boxes")
parser.add_argument("output", type=Path, help="Output directory for plots")
args = parser.parse_args()


names = {
    "beampipe": "Beam pipe",
    "sstrips": "Short Strips",
    "lstrips": "Long Strips",
    "pixel": "Pixel",
    "solenoid": "Solenoid",
    "ecalbarrel": "EMCal barrel",
    "ecalendcap": "EMCal endcap",
}


colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

with args.config.open() as f:
    config = json.load(f)


zmax = float("-inf")
zmin = float("inf")
rmax = 0

fig, ax = plt.subplots()
ax.set_xlabel("z [mm]")
ax.set_ylabel("r [mm]")

labels = set()

for color, (key, boxes) in zip(colors, config.items()):

    for box in boxes:
        label = names[key] if key not in labels else None
        labels.add(key)
        ax.add_patch(Rectangle((box["zmin"], box["rmin"]), 
                               box["zmax"]-box["zmin"], 
                               box["rmax"]-box["rmin"], fc=color, label=label))

        zmax = max(zmax, box["zmax"])
        zmin = min(zmin, box["zmin"])
        rmax = max(rmax, box["rmax"])

    print(key)

ax.set_xlim(zmin, zmax)
ax.set_ylim(0, 1.25*rmax)
ax.legend(ncol=3)


fig.tight_layout()
fig.savefig(str(args.output))

