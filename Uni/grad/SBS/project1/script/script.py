import requests
import pandas as pd
import urllib
import urllib3
import geopy.distance
from numpy import sqrt, log10, pi, exp, tanh
from math import ceil
from time import sleep

## Dictionary of possible transmitters: [x, y, max_antenna_height]
# max_antenna_height = 10 for unknown values
TRANSMITTERS = {
    "Praha-Slivenec": [50.02011399999999, 14.34015700000002, 10],
    "Bubovice": [49.9724979999999, 14.17755099999998, 10],
    "Zlicin": [50.05122599999999, 14.28301700000001, 10],
    "Rudna": [50.0186111, 14.1958333, 10],
    "Chyne": [50.0523211, 14.2266106, 10],
    "Strahov": [50.0797222, 14.3758333, 60],
    "Zahrabska": [49.9711111, 14.1116667, 10],
    "FEL": [50.1024983, 14.3927758, 40],
    "Beroun": [49.9627778, 14.0650000, 14],
}

### Simple persistent request function
def get_elevation(point):
    url = 'https://api.opentopodata.org/v1/eudem25m?'
    query = {'locations': f"{point[0]},{point[1]}"}
    while True:
        try:
            response = requests.get((url + urllib.parse.urlencode(query)))
        except (OSError, urllib3.exceptions.ProtocolError) as ex:
            print('*' * 20, 'Error Occured', '*' * 20)
            print(ex)
            continue
        if not 'results' in response.json():
            sleep(0.5)
            continue
        break
    return response.json()['results'][0]['elevation']

### Function for calculating and printing path-clearance parameters
def path_clearance(transmitter1, transmitter2, Txh = 0):
    ## Initialization
    start = TRANSMITTERS[transmitter1]
    stop =  TRANSMITTERS[transmitter2]
    Rx = transmitter2
    if get_elevation(start[0:2]) > get_elevation(stop[0:2]): start, stop, Rx = stop, start, transmitter1
    if Txh == 0: Txh = start[2]

    ## Inspect the path of the radio link for elevation every 100 meters
    d = geopy.distance.geodesic(start[0:2], stop[0:2])
    number_of_points = ceil(d.m/100)
    # number_of_points = 20     # small number of points for debugging

    lat = []
    lon = []
    difference = [stop[0] - start[0], stop[1] - start[1]]
    for i in range(number_of_points + 1):
        lat.append(start[0] + i*difference[0]/number_of_points)
        lon.append(start[1] + i*difference[1]/number_of_points)

    data = pd.DataFrame({'lat': lat, 'lon': lon})
    data['elevation'] = data.apply(get_elevation, axis=1)
    data_csv = data.to_csv(header=None, lineterminator='\n')

    ## Process the data to find the point of maximum terrain elevation
    data_rows = data_csv.split('\n')
    number_of_rows = len(data_rows) - 1
    elevations = [None]*number_of_rows
    for i in range(number_of_rows):
        elevations[i] = data_rows[i].split(',')[3]
    index_max = max(range(len(elevations))[1:-1], key=elevations.__getitem__)
    point_max = data_rows[index_max].split(',')[1:3]

    ## Calculate the radius of the first Fresnel zone F1
    f = 18                                                  # frequency of the radio link in GHz
    d1 = geopy.distance.geodesic(start[0:2], point_max)     # distance from start to the point of max elevation
    d2 = geopy.distance.geodesic(point_max, stop[0:2])      # distance from the point of max elevation to stop
    assert d - (d1 + d2) < 1e-04                            # adding the previous two should yield overall
    F1 = 17.3*sqrt((d1.km*d2.km)/(f*d.km))

    ## Calculate the subrefraction height correction x for normal atmosphere (k = 4/3)
    R_Z = 6371e3                                            # mean value of the Earth's radius in m
    k_e = 4/3                                               # value of the subrefractive factor k for normal atmosphere
    R_e = k_e*R_Z                                           # effective value of Earth's radius due to subrefraction
    x = (d1.m*d2.m)/(2*R_e)
    ## Obtain antenna height Rxh1 for full F1 clearance
    h1 = float(elevations[0]) + Txh                         # fixed antenna height in start including terrain elevation
    h0 = float(elevations[index_max])                       # maximum terrain elevation along the path
    Rxh1 = d.m*(h0 + x + F1 - h1)/d1.m + h1 - float(elevations[-1])
    ## Calculate average terrain height
    h_avg = 0
    for elevation in elevations: h_avg += float(elevation)
    h_avg = h_avg/len(elevations)

    ## Calculate the subrefraction height correction x for the path length in question
    k_e = 157/(144+2670/d.m)
    R_e = k_e*R_Z
    x = (d1.m*d2.m)/(2*R_e)
    ## Obtain antenna height Rxh2 for temperature climate with path obstriction extended along a portion of the path (0.3*F1)
    Rxh2 = d.m*(h0 + x + 0.3*F1 - h1)/d1.m + h1 - float(elevations[-1])

    ## Choose the larger of the antenna heights
    Rxh = max(Rxh1, Rxh2)

    ## Estimate diffraction loss (valid for losses greater than 15 dB) due to the path in question
    h = float(elevations[0]) + index_max*(float(elevations[-1]) - float(elevations[0]))/number_of_rows - h0
    Ad = -20*h/F1 + 10

    ## Free space loss (FSL)
    Gr = 44
    Gt = 44
    c = 3e8
    FSL = 20*log10(4*pi*f*1e9*d.m/c) - Gr - Gt

    ## Attenuation due to atmospheric gases
    gamma = 0.1                                             # rough value of gamma obtained from graphs in ITU-R P.676
    Aa = gamma*d.km

    ## Outages (Tonda)
    zasahFresnel = 1 # Pokud fresnel nepotka prekazku, dej 0
    A = 15.81 # [dB] fade depth teda vlastne margin co mame? 
    K = 1.625923e-5 #exp logK z r.530
    dn75 = 7.7142 # row 300 collumn 
    print(K)
    eps_p = abs(Txh-Rxh)/d.km
    hc= (Txh+Rxh)/2 - (d.km**2)/102 - sum([Txh,Rxh])/len([Txh,Rxh])
    vsrlimit = float((dn75*d.km**1.5 *f**0.5)/24730)
    vsr=float(min((dn75/50)**1.8 * exp(-hc/(2.5*sqrt(d.km))),vsrlimit))
    if ~zasahFresnel:
        vsr = 0
    pw = K*(d.km**3.51)*(f**2 +13)**0.447 * 10**(-(0.376*tanh((hc-147)/125))-0.334*eps_p**0.39- 0.00027*min(Txh,Rxh)+17.85*vsr-A/10)
    print(pw)

    ## Output the results into a file
    with open(fr"C:\\Workspace\\Uni\\grad\\SBS\\script\\{transmitter1}_{transmitter2}_{number_of_points}points.txt", 'w') as file:
        file.write('*'*20 + "PATH-CLEARANCE" + '*'*20 + "\n")
        file.write("Point of maximum elevation h0: [" + point_max[0] + "," + point_max[1] + "]\n")
        file.write("Parameters: d = " + str(round(d.km,2)) + " km, d1 = " + str(round(d1.km,2)) + " km, d2 = " + str(round(d2.km,2)) + " km, h1 = (" + str(round(float(elevations[0]))) + " + " + str(round(start[2])) + ") m, h0 = " + str(round(h0,2)) + " m, F1 = " + str(round(F1,2)) + " m, x = " + str(round(x,2)) + " m\n")
        file.write("Average terrain height: h_avg = " + str(round(h_avg, 2)) + " m.\n")
        file.write("Height of the " + Rx + " antenna required for line-of-sight path clearance: Rxh = " + str(round(Rxh,2)) + " m\n")
        file.write("Corresponding losses: Ad = " + str(round(Ad, 2)) + " dB, FSL = " + str(round(FSL, 2)) + " dB, Aa = " + str(round(Aa, 2)) + " dB\n")
        file.write("Outages from Tonda's part of the script: K = " + str(K) + ", pw = " + str(pw) + ".\n")
        file.write("\n" + '*'*20 + "DATA" + '*'*20 + '\n' + data_csv)

    return Rxh


### Main function
if __name__ == "__main__":
    ## Path clearance function with user input
    path_clearance("FEL", "Strahov")
    Rxh = path_clearance("Strahov", "Rudna")
    Rxh = path_clearance("Rudna", "Zahrabska", Rxh)
    path_clearance("Zahrabska", "Beroun", Rxh)
