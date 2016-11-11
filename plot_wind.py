

import matplotlib.pyplot as plt
from scipy import ndimage
import numpy as np
from pyproj import Proj, transform
import matplotlib.cm as cm

inProj = Proj(init='epsg:4326')
outProj = Proj(init='epsg:3857')
wfac  = 1.
wind_thres = 4.
fig_w = 10

path_to_card = '/media/topiko/64GbMCard/fmi_data'

fix_coord_dict = {'uts':[[587.96, 845.51], [69.9090, 27.0197], [3007819.24538701, 11039161.72817472]],
                'hel':[[516,48], [60.16393, 24.89984], [2771837.50963399, 8436325.88087744]],
                'osl':[[32,30], [59.8737, 10.64469], [1184961.47045225, 8371672.10795351]]}

def coords_toPxlCoords(coords):

    coord_mat = np.zeros((3,3))
    pix_vec_x = np.zeros(3)
    pix_vec_y = np.zeros(3)

    i = 0
    for fix in fix_coord_dict:
        coord_mat[i,:2] = fix_coord_dict[fix][2]
        pix_vec_x[i] = fix_coord_dict[fix][0][0]
        pix_vec_y[i] = fix_coord_dict[fix][0][1]
        i += 1
    coord_mat[:,2] = 1
    fac_x = np.linalg.solve(coord_mat, pix_vec_x)
    fac_y = np.linalg.solve(coord_mat, pix_vec_y)
    pxl_coords = np.zeros(np.shape(coords))

    pxl_coords[:,0] = coords[:,0]*fac_x[0] + coords[:,1]*fac_x[1] + fac_x[2]
    pxl_coords[:,1] = coords[:,0]*fac_y[0] + coords[:,1]*fac_y[1] + fac_y[2]

    return pxl_coords

def index_cities(cities, ncities):
    city_coords = np.zeros((ncities, 2))
    city_data = []
    i = 0
    for key in cities:
        city_data.append(cities[key][1:])
        lat, lon = cities[key][0] #cities_dict[city]
        city_coords[i] = transform(inProj, outProj, lon, lat)
        #city_coords[i] = cities[key][0]
        i += 1
    return np.array(city_data), city_coords

def plot_bestSpot(plot_cities, ncities):

    city_data, city_coords = index_cities(plot_cities, ncities)
    pxl_coords = coords_toPxlCoords(city_coords)
    percent_threshold = .8

    energies = city_data[:,0]
    act_percent = city_data[:,1]
    wdirs = city_data[:,2]
    plt.figure(figsize = (fig_w, fig_w*1.2))
    im = plt.imread('fin_back.png')
    im = np.flipud(im)
    #im = ndimage.rotate(im, 180, reshape = False)
    implot = plt.imshow(im, origin = 'lower')
    plt.scatter([587,516,32],[845,48,30], marker = 'o')
    plt.scatter(pxl_coords[:,0],pxl_coords[:,1], marker = 'o', color = 'black')


    for i in range(len(wdirs)):
        size = 1500*(energies[i] / np.max(energies))
        #print(energies)
        if act_percent[i] > percent_threshold:
            color_act = 'green'
        else:
            color_act = 'red'

        plt.scatter(pxl_coords[i,0],pxl_coords[i,1], marker = 'o',
            s = size, color = color_act, alpha = .6)


        if np.isfinite(wdirs[i]):
            phi = wdirs[i]/360.*2*np.pi
            x = np.sin(phi)
            y = np.cos(phi)
            dir_vec = 10*np.array([x,y])/np.sqrt(x**2 + y**2)
            plt.arrow(pxl_coords[i,0],pxl_coords[i,1],
                dir_vec[0], dir_vec[1], lw = 3, head_width=10,
                head_length=10, fc='k', ec='k', alpha = 1.)


    plt.gca().set_axis_off()
    plt.subplots_adjust(top = 1, bottom = 0, right = 1., left = .0,
            hspace = 0, wspace = 0)
    plt.margins(0.,0.)
    plt.gca().xaxis.set_major_locator(plt.NullLocator())
    plt.gca().yaxis.set_major_locator(plt.NullLocator())
    plt.gca().text(10, 900, 'Best spot: circle size correlates \n' \
            + 'with expected energy production.', fontsize = 15)
    plt.gca().text(10, 850, 'Color green: operating time > %i %%' \
        %int(100*percent_threshold), fontsize = 15)



    #plt.title(day[1], location = (0,0,.4,.1))
    plt.savefig(path_to_card + "/Best_spot.png", dpi = 100) #, bbox_inches = 'tight', pad_inches = 0)
    plt.show()
    plt.clf()
    plt.close('all')


def plot_wind(plot_cities, ncities, day = [0, '']):


    city_data, city_coords = index_cities(plot_cities, ncities)
    pxl_coords = coords_toPxlCoords(city_coords)

    city_winds = city_data[:,1]
    city_wind_dirs = city_data[:,2]

    plt.figure(figsize = (fig_w, fig_w*1.2))
    im = plt.imread('fin_back.png')
    im = np.flipud(im)
    #im = ndimage.rotate(im, 180, reshape = False)
    implot = plt.imshow(im, origin = 'lower')
    plt.scatter([587,516,32],[845,48,30], marker = 'o')
    plt.scatter(pxl_coords[:,0],pxl_coords[:,1], marker = 'o', color = 'black')


    for i in range(ncities):
        if np.isfinite(city_winds)[i]:
            windms = city_winds[i]
            size = 5*( windms / wind_thres )**(3./2)
            if wind_thres < windms:
                color_ar = 'green'
                #plt.scatter(pxl_coords[i,0],pxl_coords[i,1], marker = 'o',
                #        s = size, color = 'green', alpha = .6)
            else:
                color_ar = 'red'
                #plt.scatter(pxl_coords[i,0],pxl_coords[i,1], marker = 'o',
                #        s = size, color = 'red', alpha = .6)
            if np.isfinite(city_wind_dirs)[i]:
                phi = city_wind_dirs[i]/360.*2*np.pi
                x = np.sin(phi)
                y = np.cos(phi)
                dir_vec = size*np.array([x,y])/np.sqrt(x**2 + y**2)
                plt.arrow(pxl_coords[i,0],pxl_coords[i,1],
                    dir_vec[0], dir_vec[1], lw = size/2., head_width=10,
                    head_length=10, alpha = np.min([1, windms/12.]),
                    color=color_ar) #color = color_ar) #, fc='k', ec='k'

    plt.gca().set_axis_off()
    plt.subplots_adjust(top = 1, bottom = 0, right = 1., left = .0,
            hspace = 0, wspace = 0)
    plt.margins(0.,0.)
    plt.gca().xaxis.set_major_locator(plt.NullLocator())
    plt.gca().yaxis.set_major_locator(plt.NullLocator())
    if day[1] != '':
        plt.gca().text(10, 900, day[1], fontsize = 15)
        #plt.title(day[1], location = (0,0,.4,.1))
        plt.savefig(path_to_card + "/pics/wind%04d.png" %day[0], dpi = 50) #, bbox_inches = 'tight', pad_inches = 0)
    #plt.show()
    plt.clf()
    plt.close('all')
