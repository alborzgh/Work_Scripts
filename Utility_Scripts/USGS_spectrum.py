"""
# ########################################################### #
#                 Golder Associates, Inc.                     #
# ########################################################### #

This file includes functions to retrieve spectral acceleration 
values from the USGS website and generate design response 
spectrum based on AASHTO design manual.

Golder Associates, Inc.
Created by: Alborz Ghofrani
            Feb 2019
"""


# import required modules
import numpy as np
import matplotlib.pyplot as plt
import requests

from scipy.interpolate import interp2d
from scipy.interpolate import interp1d

def get_USGS_hazard(lat, lon, return_period = 2475):
    """
    Retrieves spectral recurrence values from the USGS
    deaggregation (B/C boundary) website and interpolates 
    at the return_period. This can be used to generate the 
    design response spectrum based on AASHTO manual.

    Parameters
    ----------
    lat, lon : are the latitude and longitude of the location
            of interest.
    return_period : is the return period at which spectral
            values are interpolated.

    Returns PGA, Sa at 0.2 sec period and Sa at 1.0 sec period.
    """

    req_link = "https://earthquake.usgs.gov/hazws/staticcurve/1/E2014R1/COUS0P05/{1}/{0}/any/760".format(lat,lon)
    r = requests.get(req_link)

    json_data = r.json()

    ground_motion = {'PGA':{}, 'SA0P2':{}, 'SA1P0':{}}
    if json_data['status'] == "success":
        res_date = json_data['date']
        data = json_data['response']
        for data_i in data:
            curve_name = data_i['metadata']['imt']['value']
            curve_display_name = data_i['metadata']['imt']['display']
            x_vals = data_i['metadata']['xvals'][:]

            lat1 = data_i['data'][0]['latitude']
            lon1 = data_i['data'][0]['longitude']
            y1_vals = np.array(data_i['data'][0]['yvals'])

            lat2 = data_i['data'][1]['latitude']
            lon2 = data_i['data'][1]['longitude']
            y2_vals = np.array(data_i['data'][1]['yvals'])

            lat3 = data_i['data'][2]['latitude']
            lon3 = data_i['data'][2]['longitude']
            y3_vals = np.array(data_i['data'][2]['yvals'])

            lat4 = data_i['data'][3]['latitude']
            lon4 = data_i['data'][3]['longitude']
            y4_vals = np.array(data_i['data'][3]['yvals'])

            int_x = np.array([lat1,lat2,lat3,lat4])
            int_y = np.array([lon1,lon2,lon3,lon4])

            interpolated_yvals = np.zeros(y1_vals.shape)

            
            for ii, _ in enumerate(interpolated_yvals):
                z = np.array([y1_vals[ii], y2_vals[ii], y3_vals[ii], y4_vals[ii]])
                f = interp2d(int_x, int_y, z, kind='linear')
                interpolated_yvals[ii] = f(lat, lon)

            calc_gm_rp = interp1d(np.log(interpolated_yvals), np.log(x_vals), kind='linear')
            
            ground_motion[curve_name]['xvals'] = x_vals
            ground_motion[curve_name]['yvals'] = interpolated_yvals
            ground_motion[curve_name]['display'] = curve_display_name
            ground_motion[curve_name]['Sa'] = np.exp(calc_gm_rp(np.log(1.0/return_period)))

    else:
        print("Some error occured!\n")
        return
    
    return (ground_motion['PGA']['Sa'], ground_motion['SA0P2']['Sa'], ground_motion['SA1P0']['Sa'], ground_motion, json_data)

def get_USGS_hazard_general(lat, lon, return_period = 2475, site_class='B', version = 2014):
    """
    Retrieves spectral recurrence values from the USGS
    deaggregation (chosen site class) website and interpolates 
    at the return_period. This can be used to generate the 
    design response spectrum based on AASHTO manual.

    Parameters
    ----------
    lat, lon : are the latitude and longitude of the location
            of interest.
    return_period : is the return period at which spectral
            values are interpolated.

    Returns PGA, Sa at 0.2 sec period and Sa at 1.0 sec period.
    """

    site_class_vs30 = {
        'A' : 2000,
        'AB': 1500,
        'B' : 1150,
        'BC': 760,
        'C' : 537,
        'CD': 360,
        'D' : 259,
        'DE': 180,
        'E' : 150
    }

    req_link = "https://earthquake.usgs.gov/nshmp-haz-ws/hazard/E{2}/COUS/{1}/{0}/any/{3}".format(lat, lon, version, site_class_vs30[site_class])
    r = requests.get(req_link)

    json_data = r.json()

    ground_motion = {}
    if json_data['status'] == "success":
        periods = []
        spectrum = []
        res_date = json_data['date']
        data = json_data['response']
        for data_i in data:
            curve_name = data_i['metadata']['imt']['value']
            curve_display_name = data_i['metadata']['imt']['display']
            x_vals = data_i['metadata']['xvalues'][:]

            if curve_name == 'PGA':
                periods.append(0)
            elif curve_name.startswith('SA'):
                periods.append(float(curve_name.replace('SA','').replace('P','.')))
            else:
                print(f"What is {curve_name}? Aborting.")
                exit()

            interpolated_yvals = np.array(data_i['data'][0]['yvalues'])

            calc_gm_rp = interp1d(np.log(interpolated_yvals), np.log(x_vals), kind='linear')
            
            ground_motion[curve_name] = {}
            ground_motion[curve_name]['xvals'] = x_vals
            ground_motion[curve_name]['yvals'] = interpolated_yvals
            ground_motion[curve_name]['display'] = curve_display_name
            ground_motion[curve_name]['Sa'] = np.exp(calc_gm_rp(np.log(1.0/return_period)))
            spectrum.append(ground_motion[curve_name]['Sa'])

    else:
        print("Some error occured!\n")
        return
    
    return (periods, spectrum, ground_motion, json_data)


def adjust_for_site_class(PGA, Ss, S1, site_class='C'):
    """
    Returns the site-class-corrected values based on AASHTO
    design manual.
    """

    if site_class == 'F':
        print('You should do a site-specific analysis!\n')
        return

    if not site_class in ['A','B','C','D','E']:
        print('Choose a valid site class!\n')
        return


    pga_x = np.array([0.1,0.2,0.3,0.4,0.5])
    Ss_x  = np.array([0.25,0.5,0.75,1.0,1.25])
    S1_x  = np.array([0.1,0.2,0.3,0.4,0.5])

    pga_Ss_fact = {'A':np.array([0.8,0.8,0.8,0.8,0.8]),
                   'B':np.array([1.0,1.0,1.0,1.0,1.0]),
                   'C':np.array([1.2,1.2,1.1,1.0,1.0]),
                   'D':np.array([1.6,1.4,1.2,1.1,1.0]),
                   'E':np.array([2.5,1.7,1.2,0.9,0.9]),                
                  }

    S1_fact =  {'A':np.array([0.8,0.8,0.8,0.8,0.8]),
                'B':np.array([1.0,1.0,1.0,1.0,1.0]),
                'C':np.array([1.7,1.6,1.5,1.4,1.3]),
                'D':np.array([2.4,2.0,1.8,1.6,1.5]),
                'E':np.array([3.5,3.2,2.8,2.4,2.4]),                
                }

    pga_func = interp1d(pga_x, pga_Ss_fact[site_class], kind='linear', fill_value=(pga_Ss_fact[site_class][0],pga_Ss_fact[site_class][-1]), bounds_error=False)
    Ss_func  = interp1d(Ss_x,  pga_Ss_fact[site_class], kind='linear', fill_value=(pga_Ss_fact[site_class][0],pga_Ss_fact[site_class][-1]), bounds_error=False)
    S1_func  = interp1d(S1_x,  S1_fact[site_class],     kind='linear', fill_value=(S1_fact[site_class][0],S1_fact[site_class][-1]),         bounds_error=False)

    pga_factor = pga_func(PGA)
    Ss_factor  = Ss_func(Ss)
    s1_factor  = S1_func(S1)

    return(pga_factor, Ss_factor, s1_factor)

def adjust_for_site_class_WSDOT_2019(PGA, Ss, S1, site_class='C'):
    """
    Returns the site-class-corrected values based on WSDOT GDM 2019
    design manual.
    """

    if site_class == 'F':
        print('You should do a site-specific analysis!\n')
        return

    if not site_class in ['A','B','C','D','E']:
        print('Choose a valid site class!\n')
        return


    pga_x = np.array([0.1,0.2,0.3,0.4,0.5,0.6])
    Ss_x  = np.array([0.25,0.5,0.75,1.0,1.25, 1.5])
    S1_x  = np.array([0.1,0.2,0.3,0.4,0.5, 0.6])

    pga_fact = {'A':np.array([0.8,0.8,0.8,0.8,0.8,0.8]),
                'B':np.array([0.9,0.9,0.9,0.9,0.9,0.9]),
                'C':np.array([1.3,1.2,1.2,1.2,1.2,1.2]),
                'D':np.array([1.6,1.4,1.3,1.2,1.1,1.1]),
                'E':np.array([2.4,1.9,1.6,1.4,1.2,1.1]),                
                  }

    Ss_fact =  {'A':np.array([0.8,0.8,0.8,0.8,0.8,0.8]),
                'B':np.array([0.9,0.9,0.9,0.9,0.9,0.9]),
                'C':np.array([1.3,1.3,1.2,1.2,1.2,1.2]),
                'D':np.array([1.6,1.4,1.2,1.1,1.0,1.0]),
                'E':np.array([2.4,1.7,1.3,1.0,0.9,0.9]),                
                }

    S1_fact =  {'A':np.array([0.8,0.8,0.8,0.8,0.8,0.8]),
                'B':np.array([0.8,0.8,0.8,0.8,0.8,0.8]),
                'C':np.array([1.5,1.5,1.5,1.5,1.5,1.4]),
                'D':np.array([2.4,2.2,2.0,1.9,1.8,1.7]),
                'E':np.array([4.2,3.3,2.8,2.4,2.2,2.0]),                
                }

    pga_func = interp1d(pga_x, pga_fact[site_class], kind='linear', fill_value=(pga_fact[site_class][0],pga_fact[site_class][-1]), bounds_error=False)
    Ss_func  = interp1d(Ss_x,  Ss_fact[site_class],  kind='linear', fill_value=(Ss_fact[site_class][0],Ss_fact[site_class][-1]),   bounds_error=False)
    S1_func  = interp1d(S1_x,  S1_fact[site_class],  kind='linear', fill_value=(S1_fact[site_class][0],S1_fact[site_class][-1]),   bounds_error=False)

    pga_factor = pga_func(PGA)
    Ss_factor  = Ss_func(Ss)
    s1_factor  = S1_func(S1)

    return(pga_factor, Ss_factor, s1_factor)




def main1():
    """ Example use of functions in this module."""


    # site location
    lat = 47.572260
    lon = -122.347509

    return_period = 150.0
    PGA, Ss, S1, _ = get_USGS_hazard(lat, lon, return_period=return_period)
    pga_factor, Ss_factor, s1_factor = adjust_for_site_class(PGA, Ss, S1, site_class='E')

    As = pga_factor * PGA
    SDS = Ss_factor * Ss
    SD1 = s1_factor * S1

    Ts = SD1 / SDS
    To = 0.2 * Ts

    xvals = np.array([0.0, To, 0.2, Ts])
    xvals = np.append(xvals, np.linspace(Ts, 1.0, 6))
    xvals = np.append(xvals, np.linspace(1.5, 10.0, 18))

    yvals = np.array([As,SDS,SDS,SDS])
    yvals = np.append(yvals, SD1 / xvals[4:])

    plt.figure()
    plt.plot(xvals, yvals)
    plt.xlabel("Period (sec)")
    plt.ylabel("$S_a$")
    plt.grid(which='both')
    plt.xlim(xmin=0.0)
    plt.ylim(ymin=0.0)
    plt.show()

def main2():
    # site location
    lat = 47.572260
    lon = -122.347509

    return_period = 2475.0
    B_periods, B_spectrum, _, _ = get_USGS_hazard_general(lat, lon, return_period = return_period, site_class='B', version = 2014)
    BC_periods, BC_spectrum, _, _ = get_USGS_hazard_general(lat, lon, return_period = return_period, site_class='BC', version = 2014)

    plt.figure()
    plt.plot(B_periods,  B_spectrum,  'b-', label='USGS2014 - Site Class B')
    plt.plot(BC_periods, BC_spectrum, 'r-', label='USGS2014 - B/C Boundary')
    plt.xlabel("Period (sec)")
    plt.ylabel("$S_a$")
    plt.legend()
    plt.grid(which='both', lw=0.5, linestyle=':')
    plt.xlim(xmin=0.0)
    plt.ylim(ymin=0.0)
    plt.show()

if __name__ == '__main__':
    # main1()
    main2()