from owslib.wfs import WebFeatureService
from xml.etree import ElementTree as ET
import numpy as np
import datetime
from itertools import islice
import pickle
import calendar
import time
import os, sys

# store data here:
path_to_card = '/media/topiko/64GbMCard/fmi_data'

# Connect to database:
wfs = WebFeatureService(url='http://data.fmi.fi/fmi-apikey/b0e29f6f-af2f-450f-942a-7c48c89c6b47/wfs',
                        version='2.0.0', timeout = 90)
print(wfs.identification.title)

# check possible queries from
# http://en.ilmatieteenlaitos.fi/open-data-manual-fmi-wfs-services

# see explanations of parameters from:
# http://data.fmi.fi/fmi-apikey/b0e29f6f-af2f-450f-942a-7c48c89c6b47/meta?observableProperty=observation&amp;param=n_man&amp;language=eng
# possibilities are: t2m, ws_10min, wg_10min, wd_10min, rh, td, r_1h, ri_10min, snow_aws, p_sea, vis, n_man, wawa
parameters = 't2m,p_sea,ws_10min,wg_10min,wd_10min,snow_aws'

def next_month(month, year):
    if month == 12:
        return '01/01/%s' %str(year + 1)
    else:
        return '01/%02d/%s' %(month + 1,str(year))

def nth(iterable, n, default=None):
    "Returns the nth item or a default value"
    return next(islice(iterable, n, None), default)


def get_npArr(iterable, month = 1, index = 0):
    out = open('arr_tmp.txt', 'wb')
    text = nth(iterable, 0).text
    out.write(bytes(text, 'UTF-8'))
    out.close()
    return np.loadtxt('arr_tmp.txt')

def store_dataForCity(city, data, month, year):
    npDataPath = path_to_card + '/%s/%s_%02d-%i.npy' %(city, city, month, year)
    try:
        pre_data = np.load(npDataPath)
        if pre_data[-1,2] < data[-1,2]:
            new_data = np.concatenate((pre_data, data))
            np.save(npDataPath, new_data)
        else:
            print('The data for time %s seems to exist already.' \
                %datetime.datetime.fromtimestamp(int(data[-1,2])).strftime('%d/%m/%YT%H:%M'))
    except IOError as ioe:
        print('No data for %s, create new .npy file.' %city)
        if not os.path.exists(path_to_card + '/%s' %city):
            os.makedirs(path_to_card + '/%s' %city)
        np.save(npDataPath, data)

def saveData(data, month, year):
    try:
        with open(path_to_card + '/cities.dict', 'rb') as handle:
            cities_dict = pickle.loads(handle.read())
            #print(cities_dict)
        exists = False
        ncities = 0
        for city in cities_dict:
            ncities += 1
            if all(cities_dict[city] == data[0][0:2]):
                exists = True
                break
        if not exists:
            print('New city founded! city%03d' %ncities)
            cities_dict['city%03d' %ncities] = data[0][0:2]
            with open(path_to_card + '/cities.dict', 'wb') as handle:
                pickle.dump(cities_dict, handle)

        store_dataForCity(city, data, month, year)


    except IOError as ioe:
        print('No cities dict exist create one: %s' %ioe)
        cities_dict = {'city000':data[0][0:2]}
        with open(path_to_card + '/cities.dict', 'wb') as handle:
            pickle.dump(cities_dict, handle)

def fechData(month, year):

    sstart = "01/%02d/%s" %(month, year)
    send = next_month(month, year)

    StartTime = calendar.timegm(time.strptime(sstart, '%d/%m/%Y'))
    EndTime = calendar.timegm(time.strptime(send, '%d/%m/%Y'))

    period = 60 # 120 mins
    timestep = 10 # mins
    N = int(period / timestep)

    if (EndTime - StartTime)%(period*60) == 0:
        irange = range(int((EndTime - StartTime)/(period*60)))
    else:
        raise

    #for i in irange:
    i = 0
    while i <= irange[-1]:
        starttime = StartTime + period*60*i #'2010-01-01T00:00:00Z' UNIX Time
        endtime = starttime + period*60 - timestep*60  #'2010-10-30T01:00:00Z'

        # check if the files already exist
        files_exist = [os.path.exists(path_to_card + '/.tmp/loc_data_%04d.npy' %i),
                       os.path.exists(path_to_card + '/.tmp/weather_data_%04d.npy' %i)]
        if all(files_exist):
            i += 1
            continue

        # Fech the data - this may fail from time to time.
        try:
            response = wfs.getfeature(storedQueryID='fmi::observations::weather::multipointcoverage',
                                    storedQueryParams={'bbox':'18,58,32,72,epsg:4326', # lat long degrees.. #'20,60,32,71', # 328000,7138000,370000,7544000',
                                    #'parameters':'windspeedms,WindDirection,,temperature',
                                    'parameters':parameters,
                                    #wind speed 10min average, wind_gust 10min, wdirection 10min, rain intensity 10 min mm/h, snow depth cm, cloud amount
                                    'crs':'EPSG::4326',
                                    'starttime':starttime, # '2010-01-01T00:00:00Z'
                                    'endtime':endtime, # '2010-01-01T01:00:00Z'
                                    'timestep':str(timestep)})
        except Exception as e:
            # If it fails do not fall into despair - for there is always another cahnge!
            print(e)
            continue

        out = open('data_tmp.gml', 'wb')
        out.write(bytes(response.read(), 'UTF-8'))
        out.close()

        # We gather the wanted pieces from the .gml file
        tree = ET.parse('data_tmp.gml')

        location_iterator = tree.iter('{http://www.opengis.net/gmlcov/1.0}positions')
        weather_iterator = tree.iter('{http://www.opengis.net/gml/3.2}doubleOrNilReasonTupleList')

        npArrWeather = get_npArr(weather_iterator, month)
        npArrLocation = get_npArr(location_iterator, month)

        np.save(path_to_card + '/.tmp/loc_data_%04d.npy' %i, npArrLocation)
        np.save(path_to_card + '/.tmp/weather_data_%04d.npy' %i, npArrWeather)

        hed = 'times: %s to %s, index: %i/%i' %(datetime.datetime.fromtimestamp(
                                    int(starttime)).strftime('%Y-%m-%dT%H:%M:%S'),
                                  datetime.datetime.fromtimestamp(
                                    int(endtime)).strftime('%Y-%m-%dT%H:%M:%S'),
                                    i, irange[-1])
        print(hed)
        i += 1

    # handle data into city arrays
    for i in irange:
        npArrLocation = np.load(path_to_card + '/.tmp/loc_data_%04d.npy' %i)
        npArrWeather = np.load(path_to_card + '/.tmp/weather_data_%04d.npy' %i)
        nelem = npArrLocation.shape[0]
        k = 0
        while k < nelem:
            coord = npArrLocation[k,:2]
            tmp_list = []
            for j in range(N + 1):
                if k + j < nelem:
                    if all(npArrLocation[k + j,:2] == coord):
                        tmp_vec = []
                        [tmp_vec.append(val) for val in npArrLocation[k + j]]
                        [tmp_vec.append(val) for val in npArrWeather[k + j]]
                    else:
                        break
                    tmp_list.append(tmp_vec)

            npAllData = np.array(tmp_list)
            saveData(npAllData, month, year)
            k += j

    # delete the temp files
    for the_file in os.listdir(path_to_card + '/.tmp/'):
        file_path = os.path.join(path_to_card + '/.tmp/', the_file)
        try:
            if os.path.isfile(file_path):
                os.unlink(file_path)
            #elif os.path.isdir(file_path): shutil.rmtree(file_path)
        except Exception as e:
            print(e)
    print('ready')


# Fech all data from year 2010

years = [2016]
for year in years:
    for month in range(1,11):
        fechData(int(month), year)
