import numpy as np
import matplotlib.pylab as plt

from io import BytesIO
from io import StringIO
from zipfile import ZipFile

from scipy.interpolate import interp2d
from scipy.interpolate import interp1d


site_classes = {
    'A/B':{'Vs30':1500, 'name':'AB'},
    'AB':{'Vs30':1500, 'name':'AB'},
    'B'  :{'Vs30':1080, 'name':'B'},
    'B/C':{'Vs30':760 , 'name':'BC'},
    'BC':{'Vs30':760 , 'name':'BC'},
    'C'  :{'Vs30':530 , 'name':'C'},
    'C/D':{'Vs30':365 , 'name':'CD'},
    'CD':{'Vs30':365 , 'name':'CD'},
    'D'  :{'Vs30':260 , 'name':'D'},
    'D/E':{'Vs30':185 , 'name':'DE'},
    'DE':{'Vs30':185 , 'name':'DE'},
    'E'  :{'Vs30':150 , 'name':'E'},
}

main_zip_file_address = r'C:\AlborzFiles\MyDesktop\Literature\USGS-Hazard-Map\0p01 Degree WUS Basin Map Data.zip'

def _get_hazard_curve(site_class='B', ordinate='PGA'):
    with ZipFile(main_zip_file_address, 'r') as USGS_zip_file:
        lower_zip_name = fr'0p01 Degree WUS Basin Map Data/2018_nshm_{site_classes[site_class]["name"]}_vs30_{str(site_classes[site_class]["Vs30"])}_0p01_degree_seattle_basin_maps.zip'
        lower_zip = BytesIO(USGS_zip_file.read(lower_zip_name))
        with ZipFile(lower_zip) as lower_zip_file:
            csv_address = fr'2018_nshm_{site_classes[site_class]["name"]}_vs30_{str(site_classes[site_class]["Vs30"])}_0p01_degree_seattle_basin_maps/{ordinate}/curves.csv'
            with lower_zip_file.open(csv_address, 'r') as curve_file:
                top_row = curve_file.readline().decode('utf-8').rstrip().split(',')[3:]
                hazard_x = np.array([float(x) for x in top_row])
                phantom_file = StringIO(curve_file.read().decode('utf-8'))
                data = np.loadtxt(phantom_file, delimiter=',', usecols=tuple(range(1,23)))
                lon = data[:,0]
                lat = data[:,1]
                hazard_y = data[:,2:]
                del data
                return (lat, lon, hazard_x, hazard_y)

def get_USGS_hazard_2018(lat, lon, site_class='B', return_period=2475):
    x_vals = np.array([0.0,0.01,0.02,0.03,0.05,0.075,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.75,1.0,1.5,2.0,3.0,4.0,5.0,7.5,10.0])
    y_vals = np.zeros(x_vals.shape)

    for ii, x in enumerate(x_vals):
        ordinate_text = ''
        if x == 0:
            ordinate_text = 'PGA'
        else:
            ordinate_text = 'SA' + str(x).replace('.','P')
            
        lat_list, lon_list, hazard_x, hazard_y = _get_hazard_curve(site_class, ordinate_text)
        

        loc_to_del = np.where(np.abs(lat_list - lat) > 0.02)
        lat_list = np.delete(lat_list,loc_to_del)
        lon_list = np.delete(lon_list,loc_to_del)
        hazard_y = np.delete(hazard_y,loc_to_del, 0)

        loc_to_del = np.where(np.abs(lon_list - lon) > 0.02)
        lat_list = np.delete(lat_list,loc_to_del)
        lon_list = np.delete(lon_list,loc_to_del)
        hazard_y = np.delete(hazard_y,loc_to_del, 0)

        cur_loc_hazard = np.zeros(hazard_x.shape)
        for jj, _ in enumerate(hazard_x):
            z = hazard_y[:,jj]
            f = interp2d(lat_list, lon_list, z, kind='linear')
            cur_loc_hazard[jj] = f(lat, lon)

        y_vals[ii] = np.exp(interp1d(np.log(cur_loc_hazard), np.log(hazard_x), kind='linear')(np.log(1.0/return_period)))
    return x_vals, y_vals



def main():
    print('150 year return period:')
    x_vals, y_vals = get_USGS_hazard_2018(lat=47.572260,lon = -122.347509, site_class='E', return_period=150)
    for x, y in zip(x_vals, y_vals):
        print(f'{x} {y}')

    print('2500 year return period:')
    x_vals, y_vals = get_USGS_hazard_2018(lat=47.572260,lon = -122.347509, site_class='E', return_period=2500)
    for x, y in zip(x_vals, y_vals):
        print(f'{x} {y}')


if __name__ == "__main__":
    main()