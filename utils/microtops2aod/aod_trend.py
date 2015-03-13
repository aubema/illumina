#!/usr/bin/env python
#usage: aod_trend.py aod_values.csv requested_times.csv
#    where: aod_values.csv contains timestamp and aod
#           requested_times.csv contains spectra name,
#                   start timestamp and end timestamp
#IMPORTANT aod_values.csv header must be date,aod
#IMPORTANT requested_times.csv header must be spectrum,start,end

import datetime
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.dates as mdates
from numpy import *

import logging
import sys

#dates_tmp
#aods_tmp


def aggregate_aods(dates, aods):
    aggregation_time = datetime.timedelta(hours=1)
    aggregated_aods = []
    aggregated_times = []
    interval_time = dates[0]
    interval_aods = []
    for idx in range(len(dates)):
        if dates[idx] > interval_time + aggregation_time:
            #Before starting new interval, we save last interval value
            aggregated_times.append(interval_time)
            aggregated_aods.append(min(interval_aods))

            #Starting new interval
            interval_aods = []
            interval_time = dates[idx]
            interval_aods.append(aods[idx])

        else:
            interval_aods.append(aods[idx])

    #Special case for writing last aggregated data
    aggregated_times.append(interval_time)
    aggregated_aods.append(mean(interval_aods))

    return (aggregated_times, aggregated_aods)
    

def plot_aggregated_aods(dates, aods, requested_dates=None, requested_aods=None):
    (agg_dates, agg_aods) = aggregate_aods(dates, aods)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(dates, aods, '+')
    ax.plot(agg_dates, agg_aods, 'ro')
    ax.plot(agg_dates, agg_aods, 'green')
    if requested_dates is not None and requested_aods is not None:
        ax.plot(requested_dates, requested_aods, 'o', color='green')
    ax.xaxis.set_major_locator(mdates.DayLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
    fig.savefig('aod-time2.png')
    plt.close()

def aod_at(requested_date, dates, aods):
    """Use this function with aggregated dates and aods"""
    dates_array = array(dates)
    idx_after = dates_array.searchsorted(requested_date)
    idx_before = idx_after - 1

    if idx_after == 0:
        logging.warning('Requested_date before aod measurment, ' + \
                'we will use first measurment. Use with care!')
        logging.warning('Requeste date is' + str(requested_date))
        return aods[0]

    after_before_delta = dates[idx_after] - dates[idx_before]
    
    slope = (aods[idx_after] - aods[idx_before]) / after_before_delta.seconds
    intercept = aods[idx_before]
    
    requested_before_delta = requested_date - dates[idx_before]
    
    requested_aod = intercept + slope * requested_before_delta.seconds

    return requested_aod



if __name__ == '__main__':
    aods_filename = sys.argv[1]
    requested_dates_filename = sys.argv[2]

    aods_data = mlab.csv2rec(aods_filename)
    dates = aods_data['date']
    aods = aods_data['aod']

    requested_dates_input = mlab.csv2rec(requested_dates_filename)

    (agg_dates, agg_aods) = aggregate_aods(dates, aods)

    requested_dates = []
    for prop in requested_dates_input:
        midtime = prop[1] + (prop[2] - prop[1])/2
        print str(midtime)
        requested_dates.append(midtime)
    requested_aods = []
    for rd in requested_dates:
        requested_aods.append(aod_at(rd, agg_dates, agg_aods))
    plot_aggregated_aods(dates, aods, requested_dates, requested_aods)

    requested_aods_filename = "requested_aods.csv"
    requested_aods_file = file(requested_aods_filename,'w')
    requested_aods_file.write("date,aod\n")
    for i in range(len(requested_aods)):
        requested_aods_file.write("%s,%f\n" % (requested_dates[i],
                requested_aods[i]))
    
