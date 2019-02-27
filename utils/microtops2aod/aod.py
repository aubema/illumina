#!/usr/bin/env python
#usage aod.py aod
#where aod is the requested wavelength for aod in nanometer

import sys
import logging
import math
import pylab
import matplotlib
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from numpy import *
from datetime import *

def aod_at_wavelength(wv):
    aod = aod_data.aot500 * math.pow( (500.0/wv), angstrom)
    logging.info( 'AOD at %fnm is %f' % (wv, aod) )
    logging.info('')

    return aod


if __name__ == '__main__':
    requested_wv = float(sys.argv[1])
    logging.basicConfig(level=logging.INFO)
    cimel_file = 'microtopsdata.csv'
    cimel = mlab.csv2rec(cimel_file, skiprows=2)
    wavelength = [380.0,500.0,870.0]
    
    requested_aods = []
    requested_aods_dates = []

    for tupleid in range(19,127): 
        aod_data = cimel[tupleid].view(recarray)
        aod_date = aod_data.date
        aod_time = aod_data.time.time() 
        aod_datetime = datetime(aod_date.year, aod_date.month, aod_date.day, 
                aod_time.hour, aod_time.minute, aod_time.second)

        #Timezone correction
        aod_datetime -= timedelta(hours=4)
        aod=[aod_data.aot380, aod_data.aot500, aod_data.aot870]
        signal_std = [aod_data.std380, aod_data.std500, aod_data.std870]

        if isnan(aod).any() or isnan(signal_std).any():
            logging.info('tuple %03d containt nan data' % (tupleid,))
        else:

            logging.info('Tuple %d for %s' % \
                    (tupleid,str(aod_datetime) ) )
            angstrom = -1 * (math.log(aod[1]/aod[2])) / \
                   (math.log(wavelength[1]/wavelength[2]))
            logging.info('AOD at 500nm is %f' % (aod_data.aot500,) )
            logging.info('Angstrom coefficient is %f' % angstrom)


            computed_aod = aod_at_wavelength(requested_wv)


            requested_aods_dates.append(aod_datetime)
            requested_aods.append(computed_aod)


            plt.plot(wavelength, aod, 'b+')
            plt.plot([requested_wv], [computed_aod,], 'ro')
            filename= str(aod_data.date) + '%03daod.png' % (tupleid,)
            plt.savefig(filename)
            plt.close()

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(requested_aods_dates, requested_aods, '+')
    ax.xaxis.set_major_locator(mdates.DayLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
    fig.savefig('aod-time.png')
    plt.close()

    out_filename= "microtops_interpolation_%.1f.csv" % (requested_wv,)
    out_file = file(out_filename, 'w')
    out_file.write("date,aod\n")
    for i in range(len(requested_aods)):
        out_file.write("%s,%f\n" % (str(requested_aods_dates[i],),
                requested_aods[i],)
                )

    out_file.close()



