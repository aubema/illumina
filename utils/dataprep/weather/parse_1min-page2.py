#!/usr/bin/env python
"""parse_1min-page2.py - Parse weather data

Usage: parse_1min-page2.py weather_file

Retreived data:
    datetimelocal : local time date and time
    utcdatetime : utc date and time
    pressure : average pressure in kPa
    relativehumdity: relative humidity

"""
import matplotlib.mlab as mlab
import numpy
import sys
import math
import time
import datetime

def get_rec_from_1min_page2(filename):
    converterd = {'stationid':str, 'tempinfo':str, 'preciptype':str,
                  'precipq':float, 'frozenprecipfreq':str,
                  'press1':float, 'press2':float, 'press3':float,
                  'tempdry':float, 'dewpt':float}
    weather_rec = mlab.csv2rec(filename, delimiter=' ', 
                               converterd=converterd,missing='M')

    return weather_rec

def get_datetime_from_weather_rec(weather_rec):
    rawdata = weather_rec['tempinfo']
    local_datetime = []
    utc_datetime = []
    for c_data in rawdata:
        c_time = time.strptime(c_data[3:15],'%Y%m%d%H%M')
        c_datetime = datetime.datetime(c_time.tm_year,
                                           c_time.tm_mon,
                                           c_time.tm_mday,
                                           c_time.tm_hour,
                                           c_time.tm_min)
        local_datetime.append(c_datetime)
        utc_datetime.append(c_datetime + datetime.timedelta(hours=-7))

    return (local_datetime, utc_datetime)
  
def get_pressure_from_weather_rec(weather_rec):
    press = []
    for c_tuple in weather_rec:
        press_hg = numpy.average([float(c_tuple['press1' ]),
                                        float(c_tuple['press2']),
                                        float(c_tuple['press3'])
                                 ])
        press_kPa = press_hg * 3.38638
        press.append(press_kPa)

    return press

def fahrenheit_to_celsius(temperature):
    return (5.0/9.0) * (temperature - 32.0)


def get_temp_from_weather_rec(weather_rec):
    dry_temp = []
    dew_point = []

    for c_tuple in weather_rec:
        dry_temp_F = float(c_tuple['tempdry'])
        dew_point_F = float(c_tuple['dewpt'])
        dry_temp_C = fahrenheit_to_celsius(dry_temp_F)
        dew_point_C = fahrenheit_to_celsius(dew_point_F)
        dry_temp.append(dry_temp_C)
        dew_point.append(dew_point_C)
    return (dry_temp, dew_point)

def get_relative_humidity(dry_temp, dew_point):
    """
    From http://www.faqs.org/faqs/meteorology/temp-dewpoint/

    es0 = reference saturation vapor pressure (es at a certain temp,
                                           usually 0 deg C)
    = 6.11 hPa

    T0  = reference temperature (273.15 Kelvin,  Kelvin = degree C +
                                                      273.15)

    Td  = dew point temperature (Kelvin)

    T   = temperature (Kelvin)

    lv  = latent heat of vaporization of water (2.5 * 10^6 joules
                                                       per kilogram)

    Rv  = gas constant for water vapor (461.5 joules* Kelvin / kilogram)

    e  = es0 * exp( lv/Rv * (1/T0 - 1/Td))
    es = es0 * exp( lv/Rv * (1/T0 - 1/T))
    """
    es0_hPa = 6.11
    T0_K = 273.15
    lv = 2.5e6
    T_K = dry_temp + 273.15
    Td_K = dew_point + 273.15
    Rv = 461.5

    e = es0_hPa * math.exp( lv / Rv * (1/T0_K - 1/Td_K))
    es = es0_hPa * math.exp( lv / Rv * (1/T0_K - 1/T_K))

    RH = e/es

    return RH



def get_relative_humidity_from_temps(dry_temperatures, dew_points):
    relative_humidities = []
    for i in range(len(dry_temperatures)):
        rh = get_relative_humidity(dry_temperatures[i], dew_points[i])
        relative_humidities.append(rh)

    return relative_humidities


if __name__ == '__main__':
    filename = sys.argv[1]

    weather_rec = get_rec_from_1min_page2(filename)

    pressure = get_pressure_from_weather_rec(weather_rec)

    (dry_temperature, dew_point) = get_temp_from_weather_rec(weather_rec)

    relative_humidity = get_relative_humidity_from_temps(dry_temperature,
                                                         dew_point)
    (local_datetime, utc_datetime)  = \
            get_datetime_from_weather_rec(weather_rec)

    out_file = file('weather.csv','w')
    out_file.write('localdatetime,utcdatetime,pressure,RH\n')
    for i in range(len(local_datetime)):
        out_file.write("%s,%s,%.4f,%.4f\n" % \
                       (local_datetime[i],utc_datetime[i],
                                pressure[i],relative_humidity[i]))
    out_file.close()


