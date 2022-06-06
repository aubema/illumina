#!/usr/bin/env python3

import os
from collections import OrderedDict

import illum
import numpy as np
import yaml

mie_path = os.path.dirname(illum.__path__[0]) + "Aerosol_optics/"
aerosol_types = OrderedDict(
    (
        ("inso", "insoluble aerosol"),
        ("waso", "water soluble aerosol"),
        ("soot", "soot"),
        ("ssam", "sea salt (acc)"),
        ("sscm", "sea salt (coa)"),
        ("minm", "mineral (suc)"),
        ("miam", "mineral (acc)"),
        ("micm", "mineral (coa)"),
        ("mitr", "mineral transported"),
        ("suso", "sulfate droplets"),
        ("fogr", "fog"),
    )
)
N_types = len(aerosol_types)
wavelengths = [500, 542, 614]

layers = ["CC", "CA", "CP", "U", "D", "MC", "MP", "MT", "ART", "ANT"]
RH_array = [0, 50, 70, 80, 90, 95, 98, 99]

rows = []
for layer in layers:
    for RH0 in RH_array:
        for wl in wavelengths:
            N = np.zeros(N_types)
            with open(mie_path + "combination_types") as f:
                combination_types = yaml.safe_load(f)

            index = list(aerosol_types.keys())
            for k, v in combination_types[layer].items():
                N[index.index(k)] = v

            ext_type = np.zeros(N_types)
            for j, type in enumerate(aerosol_types):
                RH = RH0 if type in ["waso", "ssam", "sscm", "suso"] else 0
                OPAC_filename = f"{mie_path}/OPAC_data/{type}{RH:02}"

                wl_prop = np.genfromtxt(
                    OPAC_filename,
                    skip_header=17,
                    skip_footer=120,
                    comments=None,
                )[:, 1:]
                pf_array = np.loadtxt(OPAC_filename)
                scat_angle = pf_array[:, 0]
                pf_array = pf_array[:, 1:].T

                ext_type[j] = np.interp(
                    wl / 1000, wl_prop[:, 0], wl_prop[:, 1]
                )

            # creating the variables of the specific combination of aerosols
            # single scattering albedo
            exta = np.sum(ext_type * N)

            # AOD
            Nmineral = 11  #
            h3 = 12
            z3 = 8
            h4 = 35
            extft = 0.296e-02
            extest = 0.214e-03
            extmin = 0.586e-02

            if (
                (layer == "CC")
                | (layer == "CA")
                | (layer == "CP")
                | (layer == "U")
            ):
                z1 = 8
                h1 = 2
                h2 = 2
            elif layer == "D":
                z1 = 2
                h1 = 6
                h2 = 6
            elif (layer == "MC") | (layer == "MP") | (layer == "MT"):
                z1 = 1
                h1 = 2
                h2 = 2
            elif layer == "ART":
                z1 = 99
                h1 = 2
                h2 = 2
            elif layer == "ANT":
                z1 = 8
                h1 = 10
                h2 = 10

            # aerosol layer
            # mineral/dust layer
            tau1 = exta * z1 * (1 - np.exp(-h1 / z1))
            # tau2=extmin*z2*(np.exp(-h1/z2)-np.exp(-h2/z2))
            # free troposhere
            tau3 = extft * z3 * (np.exp(-h2 / z3) - np.exp(-h3 / z3))
            # estratosphere
            tau4 = extest * (h4 - h3)

            tau = tau1 + tau3 + tau4  # +tau2

            if wl == 500:
                tau500 = tau
            else:
                if wl == 542:
                    ac542 = -np.log10(tau500 / tau) / np.log10(500 / wl)
                elif wl == 614:
                    ac614 = -np.log10(tau500 / tau) / np.log10(500 / wl)
                    ac = (ac542 + ac614) * 0.5
                    rows.append(layer + str(RH0) + ": " + str(ac) + "\n")

with open(mie_path + "../GUI_data/AC_std_values.txt", "w") as f:
    for row in rows:
        f.write(row)
