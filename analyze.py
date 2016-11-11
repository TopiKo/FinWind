import numpy as np
import pickle
import time
import calendar
from pyproj import Proj, transform
from plot_wind import plot_wind, plot_bestSpot
import datetime, calendar
import matplotlib.pyplot as plt
import os
from scipy.optimize import curve_fit, leastsq, least_squares

path_to_card = '/media/topiko/64GbMCard/fmi_data'

def get_wdirAv(wind_dirs, weights = None):

    if weights == None:
        weights = np.ones(len(wind_dirs))

    wind_dirs = weights[np.isfinite(wind_dirs)]*wind_dirs[np.isfinite(wind_dirs)]/360.*2*np.pi
    dirx, diry = np.cos(wind_dirs), np.sin(wind_dirs)
    dir_av_x, dir_av_y = np.average(dirx), np.average(diry)

    angle = np.arctan(dir_av_y/dir_av_x)/(np.pi * 2)*360
    if np.isfinite(angle):
        # the four cases of unit circle: note that north is 0 deg and east 90
        if 0 <= dir_av_x and 0 <= dir_av_y:
            wdir = 90. - angle
        elif dir_av_x < 0 and 0 <= dir_av_y:
            wdir = 360. - angle
        elif dir_av_x < 0 and dir_av_y < 0:
            wdir = 90 - angle
        elif 0 <= dir_av_x and dir_av_y < 0:
            wdir = 90 - angle
        else: wdir = np.nan
    else: wdir = np.nan

    return wdir

def get_av_withNans(data):
    return np.average(data[np.isfinite(data)])

def get_data_avsForPeriod(start_date, period, data):

    StartTime = calendar.timegm(time.strptime(start_date, '%d/%m/%YT%H:%M'))
    EndTime = StartTime + period
    istart = np.argmax(data[:,0] > StartTime)
    iend = np.argmin(data[:,0] < EndTime)
    if istart != 0 and iend != 0:
        temps = data[istart:iend, 1]
        temp_av = get_av_withNans(temps)
        winds = data[istart:iend, 3]
        wind_av = get_av_withNans(winds)
        wind_ds = data[istart:iend, 5]
        wind_dir_av = get_wdirAv(wind_ds)
        return temp_av, wind_av, wind_dir_av
    else:
        return np.nan, np.nan, np.nan


def get_cityWeatherAvs(times, city):

    period = calendar.timegm(time.strptime(times[1], '%d/%m/%YT%H:%M')) \
            - calendar.timegm(time.strptime(times[0], '%d/%m/%YT%H:%M'))

    #with open(path_to_card + '/cities.dict', 'rb') as handle:
    #        cities_dict = pickle.loads(handle.read())

    #averaged_cities = {}
    #n = 0

    #for city in cities_dict:
    city_weather_data = []
    for time_use in times:

        npDataPath = path_to_card + '/%s/combined.npy' %city
        try:
            city_data = np.load(npDataPath)
            temp_av, wind_av, wind_dir_av = get_data_avsForPeriod(time_use, period, city_data)
            time_s = calendar.timegm(time.strptime(time_use, '%d/%m/%YT%H:%M'))

            if any(np.isfinite([temp_av, wind_av, wind_dir_av])):
                city_weather_data.append([time_s, temp_av, wind_av, wind_dir_av])

        except FileNotFoundError as fnfe:
            continue

    if len(city_weather_data) != 0:
        return np.array(city_weather_data)
    else: return None
            #averaged_cities[city] = np.array(city_weather_data)
            #n += 1
            #if n > 40: break
    #return averaged_cities

def get_getBestArea(times):

    wind_threshold = 4.

    with open(path_to_card + '/cities.dict', 'rb') as handle:
        cities_dict = pickle.loads(handle.read())

    city_energies = {}
    ncities = 0
    for city in cities_dict:
        averaged_city = get_cityWeatherAvs(times, city)
        if averaged_city == None: continue

        print(city, 'ncities: %i' %ncities)
        time_s = averaged_city[:,0]
        winds = averaged_city[:,2]
        wind_dirs = averaged_city[:,3]

        wind_mask = np.isfinite(winds)
        time_wind = time_s[wind_mask]
        winds = winds[wind_mask]
        wind_dirs = wind_dirs[wind_mask]

        denergies = []
        nwinds = 0
        for i in range(len(winds) - 1):
            wind = winds[i]
            if np.isfinite(wind):
                nwinds += 1
                if wind_threshold < wind:
                     power = (wind/wind_threshold)**3
                     denergy = power*(time_wind[i+1] - time_wind[i])
                     denergies.append([denergy, wind_dirs[i]])

        if len(denergies) != 0:
            denergies = np.array(denergies)
            energy = np.average(denergies[:,0])
            wdir = get_wdirAv(denergies[:,1], denergies[:,0]/np.max(denergies[:,0]))
            active_percentage = len(denergies)/nwinds
            city_energies[city] = [cities_dict[city], energy, active_percentage, wdir]
            ncities += 1

    plot_bestSpot(city_energies, ncities)

def see_averageDevelopTemp(times):

    secInSideralYear = 365.256*3600*24
    init_time = calendar.timegm(time.strptime(str(2010), '%Y'))
    '''
    def temp_func(times, A, lam, phi, av_temp):
        print(A, lam, phi, av_temp)
        times -= init_time
        return A*np.sin(phi + times/secInSideralYear*2*np.pi) \
            + lam*times/secInSideralYear + av_temp
    '''
    with open(path_to_card + '/cities.dict', 'rb') as handle:
        cities_dict = pickle.loads(handle.read())

    fig, ax1 = plt.subplots(1,1)
    plt.show(False)
    plt.draw()
    ax1.set_xlabel('Time')
    ax1.set_ylabel('Temperature [c]')
    years = ['2010','2011','2012','2013','2014','2015','2016']
    xticks = [calendar.timegm(time.strptime(year, '%Y')) for year in years]

    plt.xticks(xticks, years)
    ax1.set_xlim([xticks[0],xticks[-1]])

    lams = []
    bounds = ([-np.inf, -np.inf, -np.pi, -np.inf], [np.inf, np.inf, np.pi, np.inf])

    #est_std, est_phase, est_mean = leastsq(optimize_func, [guess_std, guess_phase, guess_mean])[0]


    for city in cities_dict:
        averaged_city = get_cityWeatherAvs(times, city)
        if averaged_city == None: continue
        print(city)

        time_s = averaged_city[:,0]
        temps = averaged_city[:,1]
        winds = averaged_city[:,2]

        temp_mask = np.isfinite(temps)
        time_temp = time_s[temp_mask]
        temps = temps[temp_mask]

        if temps.shape[0]/len(times) < .6: continue
        print(len(temps), 12*5)


        temp_func = lambda x: x[0]*np.sin(x[1] + time_temp/secInSideralYear*2*np.pi) \
            + x[2]*(time_temp - init_time)/secInSideralYear + x[3] - temps



        x = least_squares(temp_func, [.0,0.,-np.pi/2, 0.], bounds = bounds)['x'] #[0]

        lams.append([x[3], x[2]])
        print('lam = %.3f, av_temp = %.2f c' %(x[2], x[3])) #(np.average(np.array(lams)), x[3]))
        ax1.plot(time_temp, temps, color = 'black', alpha = .4)
        plot_time = np.linspace(time_temp[0], time_temp[-1], 1000)
        ax1.plot(plot_time, x[0]*np.sin(x[1] + plot_time/secInSideralYear*2*np.pi) \
            + x[2]*(plot_time - init_time)/secInSideralYear + x[3], color = 'red', alpha = .4)

        fig.canvas.draw()

    plt.show()
    plt.clf()
    lams = np.array(lams)
    plt.scatter(lams[:,0], lams[:,1], marker = 'o')
    plt.show()

def make_video(times):

    period = calendar.timegm(time.strptime(times[1], '%d/%m/%YT%H:%M')) \
            - calendar.timegm(time.strptime(times[0], '%d/%m/%YT%H:%M'))

    with open(path_to_card + '/cities.dict', 'rb') as handle:
        cities_dict = pickle.loads(handle.read())


    plot_cities = {}
    n = 0
    for time_use in times:
        print(time_use)
        ncities = 0
        plot_cities = {}
        city_weather_data = []
        for city in cities_dict:
            city_weather_data = []
            npDataPath = path_to_card + '/%s/combined.npy' %city
            try:
                city_data = np.load(npDataPath)
                temp_av, wind_av, wind_dir_av = get_data_avsForPeriod(time_use, period, city_data)
                found_data = True
            except FileNotFoundError as fnfe:
                continue
            if found_data:
                ncities += 1
                city_weather_data = np.array(city_weather_data)
                plot_cities[city] = [cities_dict[city], temp_av, wind_av, wind_dir_av]

        plot_wind(plot_cities, ncities, [n, time_use])
        n += 1

    #subprocess.call('ffmpeg -framerate 10 -i wind%04d.png -c:v libx264 -r 30 -pix_fmt yuv420p out.mp4')


def combine_dataAndPlay(monthes):

    plot_cities = {}
    ncities = 0

    with open(path_to_card + '/cities.dict', 'rb') as handle:
        cities_dict = pickle.loads(handle.read())

    for city in cities_dict:
        city_weather_data = []
        all_data = []

        city_filesDict = {}
        for city_timeFile in os.listdir(path_to_card + '/%s' %city):
            if city_timeFile.split('_')[0] == city:
                time_str = city_timeFile.split('.')[0].split('_')[-1]
                begTime = calendar.timegm(time.strptime(time_str, '%m-%Y'))
                city_filesDict[city_timeFile] = begTime
        data_tmp = []
        for key in sorted(city_filesDict, key=city_filesDict.__getitem__):
            data_tmp.append(np.load(path_to_card + '/%s/' %city + key))
        combined_tmp = np.concatenate(data_tmp)
        a,b = combined_tmp.shape
        combined = np.zeros((a, b + 2))
        for i in range(a):
            time_str = datetime.datetime.fromtimestamp(int(combined_tmp[i,2])
                                    ).strftime('%Y/%m/%d/%H/%M')
            times = time_str.split('/')
            for j in range(5):
                combined[i,j] = float(times[j])
        combined[:,5:] = combined_tmp[:,3:]

        formater = '%04d   %02d   %02d   %02d   %02d   %.2f   %.2f   %.2f   %.2f   %.2f   %.2f'
        header = 'city location lat,long = %s,%s \n' %(cities_dict[city][0],cities_dict[city][1]) + \
                'year = yyyy / month = mm / day = dd / hour = hh / minute = m / temperature = c / ' +\
                'pressure = hPa / windSpeed = m/s / windGust = m/s / windDir = Deg / snowDepth = cm'
        np.savetxt(path_to_card + '/%s/combined.txt' %city, combined, header = header, fmt = formater)
        np.save(path_to_card + '/%s/combined.npy' %city, combined_tmp[:,2:])
        print(city)

        '''
        for month, year in monthes:
            npDataPath = path_to_card + '/%s/%s_%02d-%i.npy' %(city, city, month, year)

            try:
                city_data = np.load(npDataPath)
                all_data.append(city_data)

                temps = city_data[:,3][np.isfinite(city_data[:,3])]
                winds = city_data[:,5][np.isfinite(city_data[:,5])]
                wind_dirs = city_data[:,7] #[np.isfinite(city_data[:,7])]
                wdir = get_wdirAv(wind_dirs)

                city_weather_data.append([np.average(temps), np.average(winds), wdir])


            except FileNotFoundError as fnfe:
                continue
                #print('This file "%s_%02d-%i.npy" for %s does not exist..' \
                #    %(city, month, year, city))

        if len(city_weather_data) != 0:
            city_weather_data = np.array(city_weather_data)
            temps = city_weather_data[:,0][np.isfinite(city_weather_data[:,0])]
            winds = city_weather_data[:,1][np.isfinite(city_weather_data[:,1])]
            wind_dirs = city_weather_data[:,2][np.isfinite(city_weather_data[:,2])]

            temp_av = np.average(temps)
            wind_av = np.average(winds)
            wdir = get_wdirAv(wind_dirs)

            plot_cities[city] = [cities_dict[city], temp_av, wind_av, wdir] #wind_data

            all_data = np.concatenate(all_data)
            #formater = '%.6f   %.6f   %.11f   %.2f   %.2f   %.2f   %.2f   %.2f   %.2f'
            #np.savetxt(path_to_card + '/%s/combined.txt' %city, all_data, fmt = formater)
            np.save(path_to_card + '/%s/combined.npy' %city, all_data)

            ncities += 1

        #for time in city_data[:,2]:
        #    print(datetime.datetime.fromtimestamp(time).strftime('%Y-%m-%dT%H:%M:%S'))


    #plot_wind(plot_cities, ncities)
    '''

period = 14*24 #h
period_s = period * 3600.
times = []
start_time = calendar.timegm(time.strptime('1/1/2010', '%d/%m/%Y'))
days = 365*5 + 182

for i in range(int(days*24/period)):
    dateAndTime = datetime.datetime.fromtimestamp(
                                int(start_time + i*period_s)).strftime('%d/%m/%YT%H:%M')
    times.append(dateAndTime)

#get_getBestArea(times)
see_averageDevelopTemp(times)

#make_video(times)
exit()


'''
# This is needed only if new data is feched
list_times = []
for year in range(2010, 2017):
    for month in range(1,13):
        list_times.append([month,year])

combine_dataAndPlay(list_times)
'''
