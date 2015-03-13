#!/usr/bin/env python
"""mkcatfromdat.py - create angular fonction for illumina from dat files

Usage: mkcatfromdat.py dat_file1 dat_file2 ... dat_fileN

:author: Jean-Denis Giguere
:license: gpl
"""
import matplotlib.mlab as mlab
import sys

def normalize_dat(dat_filename):
    """Normalise dat filename and return record with normalized value"""
    names = ['value', 'angle']
    c_rec = mlab.csv2rec(dat_filename, names=names, delimiter=' ')
    total_value = sum(c_rec['value'])
    c_rec['value'] = c_rec['value'] / total_value

    return c_rec


def merge_dat(dat_rec_list):
    """Create a new dat file from a list of normalized dat file"""
    nb_rec = len(dat_rec_list)
    new_rec = dat_rec_list[0]
    new_rec['value'] = new_rec['value']/nb_rec
    for c_rec in dat_rec_list[1:]:
        new_rec['value'] += c_rec['value']/nb_rec

    return new_rec

if __name__ == '__main__':

    dat_filenames = sys.argv[1:]
    dat_recs = []
    for dat_filename in dat_filenames:
        dat_rec = normalize_dat(dat_filename)
        dat_recs.append(dat_rec)
    new_record = merge_dat(dat_recs)
    mlab.rec2csv(new_record, 'out.dat', delimiter=' ')
    
