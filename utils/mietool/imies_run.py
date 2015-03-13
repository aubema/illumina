#!/usr/bin/env python
"""imies-run.py

Prepare input for imies and make run.

"""
import mie_data

import matplotlib
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from numpy import *

import logging
import sys

import subprocess

class GenericInput(object):

    def inputs():
        raise NotImplementedError

    def run(self):
        cmd = [self.command_name,]
        execute = subprocess.Popen(cmd, stdin=subprocess.PIPE, 
                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        logging.info("Going to run %s" % (cmd,))
        logging.info("""Input file is \n%s""" % (''.join(self.inputs()),))

        (stdout_content, stderr_content) = \
                execute.communicate(''.join(self.inputs()))
        if len(stdout_content) > 0:
                logging.info(stdout_content)
        if len(stderr_content) > 0:
            for line in stderr_content:
                logging.error(stderr_content)
                raise EnvironmentError


class MkpsdInput(GenericInput):
    MONODISPERSED = 1
    GATES_GAUDIN_CHUMANN = 2
    LOG_NORMAL = 3
    GAMMA = 4
    MODIFIED_GAMMA = 5
    ROSIN_RAMMLER = 6
    INDEPENDANT_SIZE_BINS = 7

    INTERP_CST = 0
    INTERP_LINEAR = 1
    INTERP_DISCRET = 9

    def __init__(self):
        self.command_name='./mkpsd'
        self.rootname='rootname'
        self.distribution_size = MkpsdInput.LOG_NORMAL
        self.interpolation_type = MkpsdInput.INTERP_LINEAR
        self.size_parameter_increment = mie_data.size_parameter_increment
        self.psd_dist_size_param = [0, 50]
        self.mean_size = 0.418
        self.sigma = 0.35

    def inputs(self):
        lines = []
        lines.append(self.rootname + '\n')
        lines.append(str(self.distribution_size) + '\n')
        lines.append(str(self.interpolation_type) + '\n')
        lines.append(str(self.size_parameter_increment) + '\n')
        lines.append('%i %i' % tuple(self.psd_dist_size_param) + '\n')
        lines.append(str(self.mean_size) + '\n')
        lines.append(str(self.sigma) + '\n')

        return lines

        

class CmbpsInput(GenericInput):
    def __init__(self):
        self.command_name = './cmbps'
        self.rootname1 = 'rootname1'
        self.rootname2 = 'rootname2'
        self.ratio1 = mie_data.fine_number_density
        self.rootname_final = 'rootnameX'

    def inputs(self):
        lines = []
        lines.append(self.rootname1 + '\n')
        lines.append(str(self.ratio1) + '\n')
        lines.append(self.rootname2 + '\n')
        lines.append(self.rootname_final + '\n')

        return lines

class MkimiInput(GenericInput):
    SPHERE = 1
    COATED_SPHERE = 2

    PARALLEL = 1
    PERPENDICULAR = 2
    RANDOM =3

    def __init__(self):
        self.command_name = './mkimi'
        self.rootname = 'rootname'
        self.verbose = int(True)
        self.particle_type = MkimiInput.SPHERE
        self.sensor_wavelength = 0.584 # In micron
        self.min_angle = 0.0
        self.max_angle = 180.0
        self.angle_increment = 1.0
        self.polarization_state = MkimiInput.RANDOM
        self.refraction_index = (1.423, -3.35e-3)

    def inputs(self):
        lines = []
        lines.append(self.rootname + '\n')
        lines.append(str(self.verbose) + '\n')
        lines.append(str(self.particle_type) + '\n')
        lines.append(str(self.sensor_wavelength) + '\n')
        lines.append('%f %f %f \n' % (self.min_angle, self.max_angle,
            self.angle_increment))
        lines.append(str(self.polarization_state) + '\n')
        lines.append('%f %f \n' % (self.refraction_index[0], 
            -1*self.refraction_index[1]))

        return lines

class ImiesInput(GenericInput):
    def __init__(self):
        self.command_name = './imies'
        self.rootname = 'rootname'
        self.sub_bins = 5000

    def inputs(self):
        lines = []
        lines.append(self.rootname + '\n')
        lines.append(str(self.sub_bins) + '\n')

        return lines


class OpacInfo():
    def __init__(self):
        self.inso_file = 'inso.out'
        self.waso_file = 'waso.out'
        self.salt_file = 'salt.out'
        self.soot_file = 'soot.out'
        self.relative_humidities = [0,50,70,80,90,95,98,99]
        self.wavelength_nb = 19
        self.rh_section_begins = [25,56,87,118,149,180,211,242]
        self.col_names = ['wv', # Wavelength [micron]
                'extcoef', # Extinction coef. [1/km]
                'scacoef', # Scattering coef. [1/km]
                'abscoef', # Absorption coef. [1/km]
                'refracreal', # Refraction index (real part)
                'refracimag', # Refraction index (imaginary part)
                'opticaldepth', # aerosol optical depth
                ]
        self.inso_properties = self.get_type_properties(self.inso_file)
        self.waso_properties = self.get_type_properties(self.waso_file)
        self.salt_properties = self.get_type_properties(self.salt_file)
        self.soot_properties = self.get_type_properties(self.soot_file)

    def get_type_properties(self, type_filename):
        type_fh = open(type_filename)
        type_data = type_fh.readlines()
        properties_records = []
        for rh_section_begin in self.rh_section_begins:
            properties = []
            for line in range(self.wavelength_nb):
                current_line = type_data[rh_section_begin - 1 + line]
                wv = current_line[0:11]
                extcoef = current_line[10:20]
                scacoef = current_line[20:30]
                abscoef = current_line[30:40]
                refracreal = current_line[40:50]
                refracimag = current_line[50:60]
                opticaldepth = current_line[60:70]

                line_properties = [wv, extcoef, scacoef, abscoef,
                    refracreal, refracimag, opticaldepth]
                line_properties_float = \
                        [float(val) for val in line_properties]
                properties.append(tuple(line_properties_float))

            logging.debug(properties)
            properties_record = array(properties,
                    {'names': ('wv', 'extcoef', 'scacoef', 'abscoef',
                        'refracreal', 'refracimag', 'opticaldepth'),
                     'formats': (float32, float32, float32, float32,
                            float32, float32, float32) })
            logging.debug('Wavelength')
            logging.debug(properties_record['wv'])
            logging.debug('Wavelength end')
            properties_records.append(properties_record)

        return properties_records


def refraction_index(type, wavelength,humidity,opac_info):
    """
    Type supported: 'rural', 'urban', 'maritime'
    Wavelength [micron]
    Currently, humidity must be one of those supported by OPAC

    """
    rh_index = opac_info.relative_humidities.index(humidity)
    wv = opac_info.inso_properties[rh_index]['wv']
    after_idx = wv.searchsorted(wavelength)
    before_idx = after_idx - 1
    wv_before = wv[before_idx]
    wv_after = wv[after_idx]

    def wavelength_interpolation(wavelength, wv, refracreal, refracimag):
        after_idx =  wv.searchsorted(wavelength)
        before_idx = after_idx -1
        wv_before = wv[before_idx]
        wv_after = wv[after_idx]

        refracreal_before = refracreal[before_idx]
        refracreal_after = refracreal[after_idx]
        refracimag_before = refracimag[before_idx]
        refracimag_after = refracimag[after_idx]

        slope_real = (refracreal_after - refracreal_before) / \
                (wv_after - wv_before)
        intercept_real = refracreal_before - slope_real * wv_before

        slope_imag = (refracimag_after - refracimag_before) / \
                (wv_after - wv_before)
        intercept_imag = refracimag_before - slope_imag * wv_before

        return (wavelength* slope_real + intercept_real,
                wavelength* slope_imag + intercept_imag)



    inso_refracreal = opac_info.inso_properties[rh_index]['refracreal']
    inso_refracimag = opac_info.inso_properties[rh_index]['refracimag']
    inso_refrac = wavelength_interpolation(wavelength, wv,
            inso_refracreal, inso_refracimag)

    waso_refracreal = opac_info.waso_properties[rh_index]['refracreal']
    waso_refracimag = opac_info.waso_properties[rh_index]['refracimag']
    waso_refrac = wavelength_interpolation(wavelength, wv,
            waso_refracreal, waso_refracimag)

    soot_refracreal = opac_info.soot_properties[rh_index]['refracreal']
    soot_refracimag = opac_info.soot_properties[rh_index]['refracimag']
    soot_refrac = wavelength_interpolation(wavelength, wv,
            soot_refracreal, soot_refracimag)

    salt_refracreal = opac_info.salt_properties[rh_index]['refracreal']
    salt_refracimag = opac_info.salt_properties[rh_index]['refracimag']
    salt_refrac = wavelength_interpolation(wavelength, wv,
            salt_refracreal, salt_refracimag)


    if type == 'rural':
        inso_ratio = 0.3
        waso_ratio = 0.7
        soot_ratio = 0.0
        salt_ratio = 0.0
    elif type == 'urban':
        inso_ratio = 0.24
        waso_ratio = 0.56
        soot_ratio = 0.20
        salt_ratio = 0.00
    elif type == 'maritime':
        inso_ratio = 0.0
        waso_ratio = 0.0
        soot_ratio = 0.0
        salt_ratio = 1.0
    else:
        raise Error('Unsupported type')

    composed_refracreal = inso_ratio * inso_refrac[0] + \
            waso_ratio * waso_refrac[0] + \
            soot_ratio * soot_refrac[0] + \
            salt_ratio * salt_refrac[0]

    composed_refracimag = inso_ratio * inso_refrac[1] + \
            waso_ratio * waso_refrac[1] + \
            soot_ratio * soot_refrac[1] + \
            salt_ratio * salt_refrac[1]

    logging.info('Wavelength: %.3f micron, humidity: %.1f %%, type %s, give ' % \
            (wavelength, humidity, type) +
            'refaction index real part (%.3e) and imaginary part (%.3e)' % \
            (composed_refracreal, composed_refracimag)
            )

    return (composed_refracreal, composed_refracimag)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    opac_info = OpacInfo()

    for rh in  range(len(opac_info.relative_humidities)):
        humidity = opac_info.relative_humidities[rh]
        logging.debug ('Relative humidity : %f' % (humidity,) )
        logging.debug('Wavelength')
        logging.debug(opac_info.waso_properties[rh]['wv'])
        logging.debug('Refraction index (real)')
        logging.debug(opac_info.waso_properties[rh]['refracreal'])
        logging.debug('Refraction index (imaginary)')
        logging.debug(opac_info.waso_properties[rh]['refracimag'])

    requested_wv_nm = [340.0, 360., 370., 380., 390., 400., 410., 420., 430., 440., 450., 460., 470., 480., 490., 500., 510., 520., 530., 540., 550., 560., 570., 580., 590., 600., 610., 620., 630., 640., 650., 660., 670., 680., 690., 700., 710., 720., 730., 740., 750., 760., 436.0, 498.0, 546.0, 569.0, 616.0]
    requested_wv = [wv/1000.0 for wv in requested_wv_nm]
    requested_type = ['rural', 'urban', 'maritime']
    # requested_type = ['maritime']
    requested_rh = opac_info.relative_humidities


    for type in requested_type:
        for rh in requested_rh:
            for wv in requested_wv:
                logging.info('*****')
                logging.info('*****')
                logging.info('Run for %s, humidity %i %% at %f micron' % \
                        (type, rh, wv))
                logging.info('*****')
                logging.info('*****')
                refrac = refraction_index(type, wv, rh, opac_info)

                psd_base_rootname = '%s_RH%02d_%.4fum' % \
                        (type, rh, wv)
                psd_fine_rootname = psd_base_rootname + '_mod1'
                psd_coarse_rootname = psd_base_rootname + '_mod2'
                mkpsdfine = MkpsdInput()
                mkpsdfine.rootname = psd_fine_rootname
                mkpsdfine.mean_size = \
                    mie_data.get_size_parameter(type, rh, wv, 'fine')[0]
                mkpsdfine.sigma = \
                    mie_data.get_size_parameter(type, rh, wv, 'fine')[1]
                mkpsdfine.run()

                mkpsdcoarse = MkpsdInput()
                mkpsdcoarse.rootname = psd_coarse_rootname
                mkpsdcoarse.mean_size = \
                    mie_data.get_size_parameter(type, rh, wv, 'coarse')[0]
                mkpsdcoarse.sigma = \
                    mie_data.get_size_parameter(type, rh, wv, 'coarse')[1]
                mkpsdcoarse.run()

                cmbps = CmbpsInput()
                cmbps.rootname1 = psd_fine_rootname
                cmbps.rootname2 = psd_coarse_rootname
                cmbps.rootname_final = psd_base_rootname
                cmbps.run()

                mkimi = MkimiInput()
                mkimi.rootname = psd_base_rootname
                mkimi.sensor_wavelength = wv
                mkimi.refraction_index = refrac
                mkimi.run()

                imies = ImiesInput()
                imies.rootname = psd_base_rootname
                imies.run()



