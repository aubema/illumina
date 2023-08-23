#!/usr/bin/env python3

from xml.etree.ElementTree import Element, SubElement, tostring

import yaml


def gambons_inputs():
    # llegim dades dels fitxers d'illumina
    print("Reading Illumina input files")
    with open("inputs_params.in") as f:
        input_params = yaml.safe_load(f)
    with open("domain_params.in") as f:
        domain_params = yaml.safe_load(f)

    def coordinate_in_dd_mm_ss(coordinate):
        if coordinate < 0:
            coordinate = coordinate * -1
            neg = -1
        else:
            neg = 1
        dd = int(coordinate)
        mm = int((coordinate - dd) * 60)
        ss = (((coordinate - dd) * 60) - mm) * 60
        return neg, str(dd), str(mm), f"{ss:.2f}"

    latneg, lat_deg, lat_min, lat_sec = coordinate_in_dd_mm_ss(
        domain_params["latitude"]
    )
    lonneg, lon_deg, lon_min, lon_sec = coordinate_in_dd_mm_ss(
        domain_params["longitude"]
    )

    # computing local time
    hour = 1 - (int(float(domain_params["longitude"]) // 15))
    if hour < 0:
        hour = 12 - hour

    print("Creating Gambons configuration file")
    # contruccio del xml

    top = Element("MapaLocalConfiguration")

    parent1 = SubElement(top, "Place")

    child1 = SubElement(parent1, "name")
    child1.text = input_params["exp_name"].split("_")[0]
    child2 = SubElement(parent1, "lonDegrees")
    child2.text = lon_deg
    child3 = SubElement(parent1, "lonMinutes")
    child3.text = lon_min
    child4 = SubElement(parent1, "lonSeconds")
    child4.text = lon_sec
    child5 = SubElement(parent1, "lonEW")
    if lonneg == 1:
        child5.text = "E"
    elif lonneg == -1:
        child5.text = "W"
    child6 = SubElement(parent1, "latDegrees")
    child6.text = lat_deg
    child7 = SubElement(parent1, "latMinutes")
    child7.text = lat_min
    child8 = SubElement(parent1, "latSeconds")
    child8.text = lat_sec
    child9 = SubElement(parent1, "latNS")
    if latneg == 1:
        child9.text = "N"
    elif latneg == -1:
        child9.text = "S"
    child10 = SubElement(parent1, "heigh")
    child10.text = str(input_params["observer_elevation"])

    parent2 = SubElement(top, "Date")

    child11 = SubElement(parent2, "year")
    child11.text = "2020"
    child12 = SubElement(parent2, "month")
    child12.text = input_params["exp_name"].split("_")[2]
    child13 = SubElement(parent2, "day")
    child13.text = input_params["exp_name"].split("_")[1]
    child14 = SubElement(parent2, "hour")
    child14.text = str(hour)
    child15 = SubElement(parent2, "minute")
    child15.text = "27"
    child16 = SubElement(parent2, "second")
    child16.text = "18.0"

    parent3 = SubElement(top, "Map")

    child21 = SubElement(parent3, "band")
    child21.text = "TESS"
    child22 = SubElement(parent3, "stars")
    child22.text = "true"
    child23 = SubElement(parent3, "zodiacalLight")
    child23.text = "true"
    child24 = SubElement(parent3, "airglow")
    child24.text = "true"
    child25 = SubElement(parent3, "moon")
    child25.text = "false"
    child26 = SubElement(parent3, "sun")
    child26.text = "false"
    child27 = SubElement(parent3, "extinction")
    child27.text = "true"
    child28 = SubElement(parent3, "depth")  # healpix
    child28.text = "6"

    parent4 = SubElement(top, "Atmosphere")

    child31 = SubElement(parent4, "isOpacModel")
    child31.text = "true"
    child32 = SubElement(parent4, "aerosolType")
    child32.text = "MACL"  # input_params['aerosol_profile']
    child33 = SubElement(parent4, "relativeHumidity")
    child33.text = str(input_params["relative_humidity"])
    child34 = SubElement(parent4, "tau0")
    child34.text = str(input_params["aerosol_optical_depth"])
    child35 = SubElement(parent4, "isTauFromModel")
    child35.text = "false"
    child36 = SubElement(parent4, "lambdaRef")
    child36.text = "500"
    # child36 = SubElement(parent4, 'angExponent')
    # child36.text = ac
    child37 = SubElement(parent4, "aeronetSite")
    child37.text = "Barcelona"
    child38 = SubElement(parent4, "aodAeronetMode")
    child38.text = "Interpolation"
    child39 = SubElement(parent4, "pfnAeronetMode")
    child39.text = "DataPre"
    child40 = SubElement(parent4, "aeronetAerosolType")
    child40.text = "COPO"
    child41 = SubElement(parent4, "aeronetRelativeHumidity")
    child41.text = "90.0"
    child42 = SubElement(parent4, "airglowFact")
    child42.text = "80.0"
    child43 = SubElement(parent4, "airglowSpectrumFile")
    child43.text = "ESO_SkyCalc_100_10.dat"
    child44 = SubElement(parent4, "terrainReflectance")
    child44.text = "0.2"
    child45 = SubElement(parent4, "isAtmosphere")
    child45.text = "true"

    parent5 = SubElement(top, "Plot")  # no afecta al csv

    child46 = SubElement(parent5, "projection")
    child46.text = "AZIMUTHAL_EQUIDISTANT"
    child47 = SubElement(parent5, "azimuthOrigin")
    child47.text = "S"
    child48 = SubElement(parent5, "reflectEW")
    child48.text = "false"
    child49 = SubElement(parent5, "grid")
    child49.text = "ALTAZIMUTH"
    child50 = SubElement(parent5, "colorScale")
    child50.text = "SQC"
    child51 = SubElement(parent5, "minValue")
    child51.text = "14.0"
    child52 = SubElement(parent5, "maxValue")
    child52.text = "24.0"
    child53 = SubElement(parent5, "isContourLines")
    child53.text = "false"
    child54 = SubElement(parent5, "whiteBackground")
    child54.text = "false"

    # print(prettify(top))

    mydata = tostring(top)
    with open("illum_conf.xml", "wb") as myfile:
        myfile.write(mydata)
