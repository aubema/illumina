import math
import os
import sys

import illum
import numpy as np
import yaml
from scipy import interpolate

# READ DATA FROM INPUT FILES
illumpath = os.path.dirname(illum.__path__[0])
mie_path = illumpath + "/Aerosol_optics/"

with open("./inputs_params.in", "r") as finputs:
    params = yaml.safe_load(finputs)
    RH = params["relative_humidity"]
    if RH == "0":
        RH = "00"
    layer = []
    layer.append(params["aerosol_profile"])
    layer.append(params["layer_type"])

with open("./Inputs/integration_limits.dat", "r") as fwavelenght:
    lines = fwavelenght.readlines()
    bins = int(lines[0])
aerosol_types = [
    "inso",
    "waso",
    "soot",
    "ssam",
    "sscm",
    "minm",
    "miam",
    "micm",
    "mitr",
    "suso",
    "fogr",
]
contlayer = 0
for combination_type in layer:
    N = np.array([0] * len(aerosol_types), dtype="float")
    for i in range(len(aerosol_types)):
        N[i] = 0
    if (
        combination_type == "Manual"
        or combination_type == "MANUAL"
        or combination_type == "manual"
    ):
        if contlayer == 0:
            print("Particle concentration of the first layer (manual)")
        elif contlayer == 1:
            print("Particle concentration of the second layer (manual)")
        N[0] = float(input("Particle concentration of insoluble aerosol:"))
        N[1] = float(input("Particle concentration of water soluble aerosol:"))
        N[2] = float(input("Particle concentration of soot:"))
        N[3] = float(input("Particle concentration of sea salt (acc):"))
        N[4] = float(input("Particle concentration of sea salt (coa):"))
        N[5] = float(input("Particle concentration of mineral (suc):"))
        N[6] = float(input("Particle concentration of mineral (acc):"))
        N[7] = float(input("Particle concentration of mineral (coa):"))
        N[8] = float(input("Particle concentration of mineral transported:"))
        N[9] = float(input("Particle concentration of sulfate droplets:"))
        N[10] = float(input("Particle concentration of fog:"))
    elif combination_type == "CC":
        N[0] = 0.15
        N[1] = 2600
    elif combination_type == "CA":
        N[0] = 0.4
        N[1] = 7000
        N[2] = 8300
    elif combination_type == "CP":
        N[0] = 0.6
        N[1] = 15700
        N[2] = 34300
    elif combination_type == "U":
        N[0] = 1.5
        N[1] = 28000
        N[2] = 130000
    elif combination_type == "D":
        N[1] = 2000
        N[5] = 269.5
        N[6] = 30.5
        N[7] = 0.142
    elif combination_type == "MC":
        N[1] = 1500
        N[3] = 20
        N[4] = 3.2 * 10 ** -3
    elif combination_type == "MP":
        N[1] = 3800
        N[3] = 20
        N[4] = 3.2 * 10 ** -3
        N[2] = 5180
    elif combination_type == "MT":
        N[1] = 590
        N[3] = 10
        N[4] = 1.3 * 10 ** -3
    elif combination_type == "ART":
        N[1] = 1300
        N[0] = 0.01
        N[3] = 1.9
        N[2] = 5300
    elif combination_type == "ANT":
        N[9] = 42.9
        N[3] = 0.047
        N[8] = 0.0053
    else:
        print(
            "ERROR. Incorret aerosol combination definition. See aerosol_guide.txt"
        )
        sys.exit()
    contlayer = contlayer + 1
    for i in range(1, bins + 1):
        wl = (float(lines[i]) + float(lines[i + 1])) / 2
        # open files from each aerosol type
        j = 0
        ext_type = np.array([0] * len(aerosol_types), dtype="float")
        scat_type = np.array([0] * len(aerosol_types), dtype="float")
        pf_type = np.arange((181 * len(aerosol_types)), dtype="float")
        pf_type = pf_type.reshape(181, len(aerosol_types))
        for type in aerosol_types:
            if (
                type == "waso"
                or type == "ssam"
                or type == "sscm"
                or type == "suso"
            ):
                RH2 = str(RH)
            else:
                RH2 = "00"
            with open(mie_path + "/OPAC_data/" + type + RH2, "r") as finso:
                rows, cols = (78 - 17, 4)
                wl_prop = np.arange(rows * cols, dtype="float")
                wl_prop = wl_prop.reshape(rows, cols)
                scat_angle = np.array([0] * (198 - 86), dtype="float")
                pf_array = np.arange(((78 - 17) * (198 - 86)), dtype="float")
                pf_array = pf_array.reshape(78 - 17, 198 - 86)
                a = 0
                for i, line in enumerate(finso):
                    if i >= 17 and i < 78:
                        temp = line.split()
                        wl_prop[a, 0] = float(temp[1])  # wavelength
                        wl_prop[a, 1] = float(temp[2])  # ext coef
                        wl_prop[a, 2] = float(temp[3])  # scat coef
                        wl_prop[a, 3] = float(temp[6])  # ssa param
                        a = a + 1
                    if i >= 86:
                        temp = line.split()
                        scat_angle[i - 86] = float(temp[0])  # scattering angle
                        temp2 = temp[1:-1]
                        contwl = 0
                        for aa in temp2:
                            pf_array[contwl][i - 86] = float(
                                aa
                            )  # pf for a wl and a scat angle
                            contwl = contwl + 1

            ext_type[j] = np.interp(wl / 1000, wl_prop[:, 0], wl_prop[:, 1])
            scat_type[j] = np.interp(wl / 1000, wl_prop[:, 0], wl_prop[:, 2])
            # cubic interpolation for a given wl and angle
            pf_function = interpolate.interp2d(
                scat_angle, wl_prop[:, 0], pf_array, kind="cubic"
            )
            for i in range(181):
                pf_type[i, j] = pf_function(
                    i, wl / 1000
                )  # pf function i=angle, j=type
            j = j + 1

        # creating the variables of the specific combination of aerosols
        # single scattering albedo
        ext_total = 0
        scat_total = 0
        for i in range(len(aerosol_types)):
            ext_total = ext_total + ext_type[i] * N[i]
            scat_total = scat_total + scat_type[i] * N[i]
        w_total = scat_total / ext_total

        # phase function (angle)
        pf_raw = np.arange((181) * 2, dtype="float")
        pf_raw = pf_raw.reshape((181), 2)
        pf_raw[:, :] = 0
        norm_factor = 0
        for i in range(181):
            pf_raw[i, 0] = i * math.pi / 180
            for j in range(len(aerosol_types)):
                pf_raw[i, 1] = pf_raw[i, 1] + pf_type[i, j] * N[j]
            norm_factor = (
                norm_factor
                + (1 / (4 * math.pi))
                * 2
                * math.pi
                * pf_raw[i, 1]
                * math.sin(pf_raw[i, 0])
                * math.pi
                / 180
            )
        with open(
            "./Inputs/" + combination_type + "_" + str(int(wl)) + ".txt", "w+"
        ) as f:
            f.write(str(w_total) + " # single scatering albedo\n")
            f.write("ScatAngle PhaseFct\n")
            pf_norm = np.arange((181) * 2, dtype="float")
            pf_norm = pf_norm.reshape((181), 2)

            for i in range(181):
                pf_norm[i, 1] = pf_raw[i, 1] / norm_factor
                pf_norm[i, 0] = pf_raw[i, 0]
                f.write(str(float(i)) + " " + str(pf_norm[i, 1]) + "\n")
