#!/usr/bin/env python3

import os
import sys
import time
import zipfile
from contextlib import suppress
from glob import glob
from subprocess import call, check_output
from threading import Thread
from urllib.request import urlopen

import click
import illum
import numpy as np
import pandas as pd
import progressbar
import requests
import yaml
from PyQt5.QtWidgets import QApplication, QDialog

from .gambons_inputs import gambons_inputs
from .light_ui import Ui_ILLUMINA


class Ui_LIGHT(Ui_ILLUMINA):
    def setupUi(self, ILLUMINA):
        super().setupUi(ILLUMINA)

        self.exp_def_button.clicked.connect(self.defining_exp)
        self.add_tech_btn.clicked.connect(self.defining_combo)
        self.clean_ls_btn.clicked.connect(self.clean_ls)
        self.aerosol_opac_button.clicked.connect(self.opac_aod)
        self.Run_simulation.clicked.connect(self.run_simulation)

        # global variables
        self.switch1 = 0
        self.perc = []
        self.tech = []
        self.ulor = []
        self.inventory_line = []
        self.atm_type = ""
        self.num_batch = 0
        self.jobs = 5 * 2  # 5 capes i 2 bandwidth
        self.illumpath = os.path.dirname(illum.__path__[0])
        self.datapath = os.path.join(self.illumpath, "GUI_data")
        self.pathparent = os.getcwd()

    def defining_exp(self):
        self.log_edit.setPlainText("Defining experiment")
        self.log_edit.repaint()

        call(["illum", "init"])

        lat = self.latitude_edit.text()
        long = self.longitude_edit.text()
        # name = self.name_edit.text()
        alt = self.altitude_edit.text()
        radius = str(70)
        if self.obstacles_checkbox.isChecked():
            fobs = str(0)
            self.dobs_edit.setText("-")
            self.hobs_edit.setText("-")
            self.dobs_edit.setReadOnly(True)
            self.hobs_edit.setReadOnly(True)
            dobs = str(10)
            hobs = str(10)
        else:
            fobs = str(0.8)
            dobs = self.dobs_edit.text()
            hobs = self.hobs_edit.text()
        hlamp = self.hlamp_edit.text()

        # creating inventory
        inv = " ".join(
            f"{p}_{t}_{u}" for p, t, u in zip(self.perc, self.tech, self.ulor)
        )
        line = "\t".join([lat, long, radius, hobs, dobs, fobs, hlamp, inv])
        self.inventory_line.append(line + "\n")

        cpus = os.cpu_count()
        if cpus is None:
            cpus = 2

        jobs_batch = int(np.ceil(self.jobs / cpus))
        self.num_batch = int(np.ceil(self.jobs / jobs_batch))

        date = self.date_edit.text()
        date_day = date.split("/")[0]
        date_month = date.split("/")[1]
        # date_year = date.split('/')[2]

        array_months = [
            "January",
            "February",
            "March",
            "April",
            "May",
            "June",
            "July",
            "August",
            "September",
            "October",
            "November",
            "December",
        ]
        name_month = array_months[int(date_month) - 1]

        # managing files
        viirs = glob(f"{self.datapath}/VIIRS_database/*{date_month}*")[0]
        with suppress(FileExistsError):
            os.makedirs("VIIRS-DNB")
        with suppress(FileExistsError):
            os.symlink(viirs, "VIIRS-DNB/" + os.path.basename(viirs))
        with suppress(FileExistsError):
            os.symlink(f"{self.datapath}/SRTM", "SRTM")
        with suppress(FileExistsError):
            os.symlink(f"{self.datapath}/hydropolys.zip", "hydropolys.zip")
        with suppress(FileExistsError):
            os.symlink(
                f"{self.datapath}/spectral_bands.dat", "spectral_bands.dat"
            )

        with open("domain_params.in", "w") as f:
            f.write(
                f"latitude: {lat}\n"
                f"longitude: {long}\n"
                "srs: auto\n"
                "scale_factor: 1.666666\n"
                "nb_pixels: 17\n"
                "nb_layers: 5\n"
                "scale_min: 750\n"
                "buffer: 10\n"
            )
        # global self.inventory_line
        if os.path.isfile("inventory.txt"):
            os.remove("inventory.txt")
        with open("inventory.txt", "w") as f:
            for line in self.inventory_line:
                f.write(line)
        self.inventory_line = []
        self.progressBar.setValue(10)
        call(["illum", "domain"])
        self.progressBar.setValue(25)
        self.log_edit.setPlainText(
            "The VIIRS montly composite image of "
            + name_month
            + " of 2020 will be used as input for the radiance"
            + " emitted by light sources"
        )
        self.log_edit.repaint()
        call(["illum", "warp"])
        self.progressBar.setValue(50)

        self.atm_type = self.comboBox_atm.currentText()
        RH = self.comboBox_rh.currentText()

        aod = self.box_aod.text()
        ac = self.box_alpha.text()
        aeh = self.box_aeh.text()

        with open(self.datapath + "/inputs_params.in", "r") as f:
            ip = yaml.safe_load(f)
        if os.path.isfile("inputs_params.in"):
            os.remove("inputs_params.in")

        ip["aerosol_profile"] = self.atm_type.split(" ")[0]
        ip["observer_elevation"] = float(alt)
        ip["relative_humidity"] = float(RH)
        ip["aerosol_optical_depth"] = float(aod)
        ip["angstrom_coefficient"] = float(ac)
        ip["aerosol_height"] = float(aeh)
        ip["exp_name"] = f"{self.name_edit.text()}_{date_day}_{date_month}"

        with open("inputs_params.in", "w") as f:
            yaml.safe_dump(ip, f)

        call(["illum", "inputs"])
        self.progressBar.setValue(80)

        pathinputs = "./Inputs"
        os.chdir(pathinputs)
        call(["illum", "batches", "-N", str(jobs_batch)])
        self.progressBar.setValue(100)

        pnorm = 100 / sum(map(float, self.perc))
        log = (
            "Experiment defined\n"
            "Simulation of the diffuse radiance in zenith direction"
            " for the V band\n"
            f"Lat: {lat} Long: {long} "
            f"Altitude relative to ground: {alt}m "
            f"Date: {date_month}/2021\n"
            f"Height of the lamps: {hlamp} "
            + (
                "No obstacles\n"
                if not self.obstacles_checkbox.isChecked()
                else f"Obstacle height: {hobs} "
                f"Distance between lamps and obstacles: {dobs}\n"
            )
            + "Inventory of light sources:\n"
            + "\n".join(
                f"{float(p)*pnorm:.4g}% {t} {u}"
                for p, t, u in zip(self.perc, self.tech, self.ulor)
            )
        )
        self.log_edit.setPlainText(log)
        self.progressBar.setValue(0)

    def run_simulation(self):
        self.log_edit.setPlainText(
            "Running the simulation.\nIt might take several minutes."
        )
        self.log_edit.repaint()

        def calling_threads(number):
            call(["bash", "batch_" + number])

        threads = []
        for i in range(self.num_batch):
            x = Thread(target=calling_threads, args=(str(i + 1)))
            x.start()
            threads.append(x)

        job_cont = 0
        self.progressBar.setValue(1)
        wavel = np.loadtxt("wav.lst")  # ['507.25','545.75','584.25','622.75']
        layer = list(range(5))
        while job_cont < self.jobs:
            for wl in wavel:
                for ly in layer:
                    time.sleep(0.5)
                    path = (
                        "./exec/elevation_angle_90/azimuth_angle_0/"
                        f"wavelength_{wl}/layer_{ly}"
                    )
                    if os.path.isfile(path + "/finished.txt"):
                        job_cont += 1
                        os.remove(path + "/finished.txt")
                        self.progressBar.setValue(
                            ((100 * job_cont) // self.jobs) - 5
                        )
        for t in threads:
            t.join()

        # GAMBONS
        mag_ref = 21.93
        radiance_ref = 2.22e-7
        os.chdir(self.pathparent)
        if self.natural_light.isChecked():
            self.progressBar.setValue(0)
            self.log_edit.setPlainText("Computing natural brightness")
            gambons_inputs()

            pathgambons = self.datapath + "/GambonsV2"
            os.chdir(pathgambons)
            output_path = os.path.abspath("./output")
            for f in os.listdir(output_path):
                os.remove(os.path.join(output_path, f))
            path = f"{self.pathparent}/illum_conf.xml"
            call(
                f"java -jar {self.datapath}/GambonsV2/Gambons.jar "
                f"-cf {path}".split()
            )
            os.chdir(self.pathparent)
            csv_file = glob(os.path.join(output_path, "*.csv"))[0]
            gambons_output = pd.read_csv(
                csv_file, names=["az", "alt", "L", "mag"], skiprows=1
            )
            radiance_V_nat = gambons_output.L[
                gambons_output.alt == gambons_output.alt.min()
            ].mean()
            mag_V_nat = mag_ref - 2.5 * np.log10(radiance_V_nat / radiance_ref)
            os.rename(
                glob(os.path.join(output_path, "*.png"))[0],
                "natural_sky.png",
            )
        self.progressBar.setValue(100)

        # extracting results
        self.log_edit.setPlainText("Simulation finished\nExtracting results")
        self.log_edit.repaint()
        lines = check_output(["illum", "extract", "-c"])
        lines = sorted(lines.decode().strip().split("\n"))
        central_wl, radiance = zip(
            *[
                [float(elem.split("_")[-1]) for elem in line.split()]
                for line in lines
            ]
        )
        wl, sens = np.loadtxt("Lights/JC_V.dat", skiprows=1).T
        mask = (wl >= 470) & (wl <= 740)
        wl_filt = wl[mask]
        sens_filt = sens[mask] / 100  # % to frac

        bins = np.loadtxt(self.datapath + "/spectral_bands.dat", delimiter=",")
        bool_array = (wl_filt >= bins[:, 0, None]) & (
            wl_filt < bins[:, 1, None]
        )
        avg_value = [np.mean(sens_filt[mask]) for mask in bool_array]
        radiance_V_art = np.sum(
            np.prod([radiance, avg_value, bins[:, 1] - bins[:, 0]], 0)
        )

        contrib_filenames = [
            glob(f"*wavelength_{wl:.1f}.hdf5")[0] for wl in central_wl
        ]
        contrib = illum.MultiScaleData.from_domain("domain.ini")
        for n, fname in enumerate(contrib_filenames):
            ds = illum.MultiScaleData.Open(fname)
            for i, layer in enumerate(ds):
                contrib[i] += layer * avg_value[n]
            os.remove(fname)
        contrib.save("contribution_map")

        mag_V_art = mag_ref - 2.5 * np.log10(radiance_V_art / radiance_ref)

        log = f"Artificial sky brightness: {mag_V_art:.2f} mag/arcsec\n"
        if self.natural_light.isChecked():
            radiance_V_tot = np.mean(radiance_V_nat) + radiance_V_art
            mag_V_tot = mag_ref - 2.5 * np.log10(radiance_V_tot / radiance_ref)
            log += f"Natural sky brightness: {mag_V_nat:.2f} mag/arcsec\n"
            log += f"Total sky brightness: {mag_V_tot:.2f} mag/arcsec\n"
        else:
            mag_V_nat = np.nan
            mag_V_tot = mag_V_art

        self.log_edit.setPlainText(log)

        # Write results file
        pnorm = 100 / sum(map(float, self.perc))
        inv_text = "".join(
            f"\n    {p}% {t} with {u}% ULOR"
            for p, t, u in zip(self.perc * pnorm, self.tech, self.ulor)
        )

        output = (
            ("Experiment name", self.name_edit.text()),
            ("Latitude", self.latitude_edit.text()),
            ("Longitude", self.longitude_edit.text()),
            ("Direction", "Zenith"),
            ("Band", "V"),
            ("Date (ddmmyyyy)", self.date_edit.text()),
            ("Inventory", inv_text),
            ("Lamps height", self.hlamp_edit.text()),
            ("Street width", self.dobs_edit.text()),
            ("Bulding's number of stories", self.hobs_edit.text(), "m"),
            ("Aerosol type", self.comboBox_atm.currentText().split(" ")[0]),
            ("Relative humidity,", self.comboBox_rh.currentText(), "%"),
            ("Artificial radiance", str(round(mag_V_art, 2)), "mag/arcsec²"),
        )
        if self.natural_light.isChecked():
            output += (
                ("Natural radiance", str(round(mag_V_nat, 2)), "mag/arcsec²"),
                ("Total radiance", str(round(mag_V_tot, 2)), "mag/arcsec²"),
            )

        with open(f"Results_experiment_{self.name_edit.text()}.txt", "w") as f:
            for elem in output:
                line = f"{elem[0]}: {elem[1]}"
                if len(elem) == 3:
                    line += f" {elem[2]}"
                f.write(line + "\n")

    def defining_combo(self):
        if self.switch1 == 0:
            self.perc = []
            self.tech = []
            self.ulor = []
        self.perc.append(self.perc_edit.text())
        self.tech.append(self.comboBox_tech.currentText())
        self.ulor.append(self.comboBox_ulor.currentText())
        self.switch1 = 1
        self.perc_edit.setText("")
        self.comboBox_tech.setCurrentIndex(0)
        self.comboBox_ulor.setCurrentIndex(0)
        message = "".join(
            f"Light fixture added: {p}% {t} with {u}% ULOR\n"
            for p, t, u in zip(self.perc, self.tech, self.ulor)
        )
        self.log_edit.setPlainText(message)

    def clean_ls(self):
        self.perc = []
        self.tech = []
        self.ulor = []
        self.log_edit.setPlainText("")

    def opac_aod(self):
        self.atm_type = self.comboBox_atm.currentText()
        RH = self.comboBox_rh.currentText()

        with open(self.datapath + "/AoD_std_values.txt") as f:
            aod = yaml.safe_load(f)
        with open(self.datapath + "/AC_std_values.txt") as f:
            ac = yaml.safe_load(f)
        atm_type = self.atm_type.split(" ")[0]
        if atm_type == "D":
            aerosol_h = 6000
        elif atm_type == "ANT":
            aerosol_h = 10000
        else:
            aerosol_h = 2000

        self.box_aod.setText(str(aod[atm_type + RH])[:5])
        self.box_alpha.setText(str(ac[atm_type + RH])[:5])
        self.box_aeh.setText(str(aerosol_h))


class AppWindow(QDialog):
    def __init__(self):
        super().__init__()
        self.ui = Ui_LIGHT()
        self.ui.setupUi(self)
        self.show()


def progress_download(url, outname=None, block_size=1024 ** 2):
    site = urlopen(url)
    meta = site.info()
    response = requests.get(url, stream=True)
    bytes = int(meta["Content-Length"])

    if outname is None:
        outname = os.path.basename(url)

    print(url)
    print(outname)

    pbar = progressbar.DataTransferBar(max_value=bytes).start()
    with open(outname, "wb") as f:
        for i, data in enumerate(response.iter_content(block_size)):
            pbar.update(block_size * i)
            f.write(data)
    pbar.finish()


def download_data():
    # Download needed data
    illumpath = os.path.dirname(illum.__path__[0])
    base_url = "http://dome.obsand.org:2080/wiki/html/"

    def ensure_downloaded(name):
        if not os.path.exists(
            os.path.join(illumpath, "GUI_data", os.path.basename(name))
        ):
            zipped = name.endswith(".zip")
            if not zipped:
                name += ".zip"
            basename = os.path.basename(name)
            path = os.path.join(illumpath, "GUI_data", basename)
            print(f"`{basename}`")
            progress_download(base_url + name, path)

            if not zipped:
                with zipfile.ZipFile(path, "r") as zipObj:
                    zipObj.extractall(path[:-4])

                os.remove(path)

    print("Downloading necessary data. It may take several minutes.")
    ensure_downloaded("illumina_light/SRTM")
    ensure_downloaded("illumina_light/GambonsV2")
    ensure_downloaded("illumina_light/VIIRS_database")
    ensure_downloaded("hydropolys.zip")

    print("Done.")


@click.command()
def light():
    """Executes the Illumina Light GUI."""
    download_data()
    app = QApplication(sys.argv)
    w = AppWindow()
    w.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    light()
