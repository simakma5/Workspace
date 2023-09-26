# importing Python libraries
#matplotlib inline
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as sst
import os


def Main():
    # INPUT
    # I - Large scale effects / Empirical model
    
    Meas_name = "Jug._Partyzanu", "Terronska"
    # Pick either 1 or 2
    Meas = 1
    Meas_name = Meas_name[Meas-1]

    try:
        os.makedirs(f'figures/{Meas_name}/1_slope')
        print("Folders created!")
    except FileExistsError:
        print("Folder already exists")


    if Meas==2:
        raw_data_rx_1 = './data/RecTrace_007.csv' # data representing equidistant measurement along a straight path
        raw_data_rx_2 = './data/RecTrace_008.csv'
    
    else:
        raw_data_rx_1 = './data/RecTrace_000.csv' # data representing equidistant measurement along a straight path
        raw_data_rx_2 = './data/RecTrace_001.csv'
    

    d_start = 2 # path starting point - antennas distance (m)
    d_end = 532 # poth end-point - antennas distance (m)
    Fkonst = 58 # a value from the power balance used below to extract the loss (dB)

    cut = [10,105,532] # vzdialenosti v ktorych budem rezat
    
    # Rx = read_data(raw_data_rx_1, d_start, d_end)


    dist, Rx = zip_data(read_data(raw_data_rx_1, d_start, d_end), read_data(raw_data_rx_2, d_end, d_start), d_start, d_end)

    dist = np.array(dist)

    # generate the distance and relative loss vectors
    # - recordings at a constant speed between d_start (m) and d_end (m) considered
    # - the loss L normalized based on the power budget
    
    d = np.array(range(len(Rx)))
    d = d_start + d * (d_end-d_start)/len(d)
    L = Rx[0] - np.array(Rx) + Fkonst
    plt.figure(3)

    #plt.plot(d,L, 'o')
    plt.semilogx(dist, L, '+', lw=0.1, label='measurement')

    plt.title(f'Path Loss Measurement')
    plt.xlabel('distance (m)')
    plt.grid()
    plt.ylabel('path loss (dB)')
    plt.show(block=False)
    plt.savefig(f"figures/{Meas_name}/Path_loss_measurement", dpi= 500)
    

    Lf = moving_average(L, dist, 30)

    #PATH LOSS Modelling

    #one_slope_fit(dist, Lf, Meas_name, plot=True)



    d_matrix, L_matrix, fit_matrix = slopes_fitting(dist, d_start, d_end, Lf, cut, Meas_name)

    if len(cut) == 3:
        
        try:
            os.makedirs(f'figures/{Meas_name}/2_slope')
            print("Folders created!")
        except FileExistsError:
            print("Folder already exists")

        plt.figure(figsize=(9, 5))
        plt.semilogx(d, L, '+', lw=0.1, label='measurement')
        # treba pridat cyklus na pridavanie plotov
        plt.semilogx(d_matrix[0], fit_matrix[1][0], lw=3, label=f'1SM: Y intercept = {fit_matrix[1][1]:.1f} dB, n = {0.1*fit_matrix[1][2]:.1f} dB') 
        plt.semilogx(d_matrix[1], fit_matrix[2][0], lw=3, label=f'1SM: Y intercept = {fit_matrix[2][1]:.1f} dB, n = {0.1*fit_matrix[2][2]:.1f} dB')       
        
        plt.xlabel('distance')
        plt.ylabel('path loss (dB)')
        plt.title('Path loss prediction vs. measurement')
        plt.legend()
        plt.grid()
        plt.savefig(f"figures/{Meas_name}/2_slope/Prediction_{len(cut)-1}_slopes", dpi= 500)
        plt.show(block=True)

        # calculate prediction error and print results
        Lp = np.concatenate((fit_matrix[1][0], fit_matrix[2][0]))
        L_cut = np.concatenate((L_matrix[0], L_matrix[1]))
        d_cut = np.concatenate((d_matrix[0], d_matrix[1]))
        
        err = Lp-L_cut
        stdev = np.std(err)
        print('-----------------------------------------------------------')
        print('-----------------------------------------------------------')
        print('Empirical Two Slope Model parameters:')
        print('standard deviation =', stdev.round(1), 'dB')
        print('mean error =', np.mean(np.abs(err)).round(1), 'dB')          
        
        
        # Large-scale Fading (Shadowing) Statistics
        # PDF
        plt.figure()
        plt.hist(err, bins=40, density=True, label='Empirical PDF')
        # comparison with Log-Normal distribution
        pdfx = np.linspace(-30,30)
        ppdf = sst.norm.pdf(pdfx, scale=stdev)
        plt.plot(pdfx, ppdf, label='Normal distribution fit')
        plt.legend()
        plt.grid()
        plt.show(block=False)
        plt.savefig(f"figures/{Meas_name}/2_slope/PDF", dpi= 500)

        plt.figure(figsize=(9, 7))
        xcdf, ycdf = ecdf(err, 50)
        plt.semilogy(xcdf, 100*ycdf, label='Empirical CDF')
        plt.semilogy(xcdf,100*sst.norm.cdf(xcdf,scale=stdev), label='Normal CDF fit')
        plt.ylim(0.1,100)
        plt.xlabel('Normalized predicted RSS (dB)')
        plt.ylabel('Percentage of locations RSS < predicted RSS (%)')
        plt.legend()
        plt.grid()
        plt.show(block=False)
        plt.savefig(f"figures/{Meas_name}/2_slope/CDF", dpi= 500)

        Pn = 20 # Normalizing the power???
        plt.figure(figsize=(9, 9))
        plt.subplot(211)
        plt.semilogx(d_cut, Pn-L_cut, '+')
        plt.plot(d_cut, Pn-Lp, lw =4, label='50 % location exceeding this' )
        plt.plot(d_cut, Pn-Lp-9.5, lw =4, label='90 % location exceeding this' )
        plt.plot(d_cut, Pn-Lp+9.5, lw =4, label='10 % location exceeding this' )
        plt.xlabel('Distance (m)')
        plt.ylabel('RSS (dBm)')
        plt.title('Predicted RSS (Distance in log and linear scale)')
        plt.legend()
        plt.grid()
        plt.subplot(212)
        plt.plot(d_cut, Pn-L_cut, '+')
        plt.plot(d_cut, Pn-Lp, lw =4, label='50 % location exceeding this' )
        plt.plot(d_cut, Pn-Lp-9.5, lw =4, label='90 % location exceeding this' )
        plt.plot(d_cut, Pn-Lp+9.5, lw =4, label='10 % location exceeding this' )
        plt.xlabel('Distance (m)')
        plt.ylabel('RSS (dBm)')
        plt.grid()
        plt.savefig(f"figures/{Meas_name}/2_slope/RSS", dpi= 500)
        plt.show(block=True)


def read_data(raw_data_I_rx, d_start, d_end):
    # Reading raw measurement data 
    with open(raw_data_I_rx) as ff:
        rRX = ff.read().split() # each line = 1601 readings from the R&S PR100 receiver screen
    if d_start > d_end: # reverse the data so the distance would be encreasing
        rRX = [x for x in reversed(rRX)]
        d_start, d_end = d_end, d_start
    print(len(rRX), 'readings')
    screens_dBm = []
    for line in rRX:
        screens_dBm.append( [float(x) for x in line.split(';')])
    # plt.figure(1)
    # plt.plot(screens_dBm[0])
    # plt.title('PR100 sample frequency sweep screen - Rx power (dBm)')
    # plt.show(block=False)
    # plot vykresli vykonove frekvencne spektrum - jedno z merani na trase. 

    # zaujima nas MAX hodnota - to je vykon sledovanej frekvencie
    # extracting RX values from raw data as a maximum value of each sweep
    Rx = [max(x) for x in screens_dBm]
    # plt.figure(2)
    # plt.plot(Rx)
    # plt.title('Measured Rx power (dBm)')
    # plt.show(block=False)

    return Rx


def zip_data(data_1, data_2, d_start, d_end): # combine two sets of data by zipping them into each other based on their corresponding d value

    dist_1 = np.array(range(len(data_1)))
    dist_2 = np.array(range(len(data_2)))

    dist_1 = d_start + dist_1 * (d_end-d_start)/len(data_1)
    dist_2 = d_start + dist_2 * (d_end-d_start)/len(data_2)
    
    #dist_1 = array.array(a_np.tobytes())
    dist_1 = dist_1.tolist()
    dist_2 = dist_2.tolist()
    
    
    dist_comb = []

    cur_dist = d_start
    Rx_comb = []

    # the ZERO step of the data ZIPIN combining
    if len(data_1) > len(data_2):       # if the first data array is longer the COMBINED data starts with data from the first array
        Rx_comb.append(data_1.pop(0))   # cut the first element from data_1 and insert it into Rx_comb
        dist_comb.append(dist_1.pop(0)) # cut the first element from dist_1 and insert it into dist_comb

        Rx_comb.append(data_2.pop(0))   # cut the first element from data_2 and insert it into Rx_comb
        dist_comb.append(dist_2.pop(0)) # cut the first element from dist_1 and insert it into dist_comb

        d_end = dist_1[-1]              # this value terminates the while cycle the biggest value is needed

    elif len(data_1) < len(data_2):
        Rx_comb.append(data_2.pop(0))       
        dist_comb.append(dist_2.pop(0))

        Rx_comb.append(data_1.pop(0))      
        dist_comb.append(dist_1.pop(0))

        d_end = dist_2[-1]   

    else:                               # both arrays have the same length -> increase the dist value of the dist_2 by one
        dist_2 = np.array(len(data_2))
        dist_2 = d_start + 1 + dist_2 * (d_end-d_start)/len(data_2)
        dist_2 = dist_2.tolist()

        d_end = dist_2[-1]

    
    # combining of the data in ZIPIN fashion
    while cur_dist < d_end:
        # get differences between current dist and the dist_1 and dist_2 lowest values
        if len(dist_1) > 0: # check if there are elements left in dist_1
            dif1 = dist_1[0] - cur_dist
        else:               # when there are no elements left increase the difference to a large value
            dif1 = 5000

        if len(dist_2) > 0:
            dif2 = dist_2[0] - cur_dist
        else:
            dif2 = 5000



        if dif1 < 0 or dif2 < 0:
            print("ERR in the zipin_data function: dif < 0 !!!")
            return [],[]

        if dif2 > dif1: # dist_1 is closer to cur_dist
            Rx_comb.append(data_1.pop(0))   # cut the first element from data_1 and insert it into Rx_comb
            cur_dist = dist_1.pop(0)        # update current distance
            dist_comb.append(cur_dist)      # cut the first element from dist_1 and insert it into dist_comb

        elif dif1 > dif2: # dist_2 is closer to cur_dist
            Rx_comb.append(data_2.pop(0))
            cur_dist = dist_2.pop(0)
            dist_comb.append(cur_dist)

        else:
            print("ERR in the zipin_data function difs are EQUAL!!!")
            return [],[]
        
    return dist_comb, Rx_comb


def moving_average(L,d, win):
    # slide window shoud consider the wavelength and environment
    Lf=[]
    for i in range(len(L)):
        imin = max(0,i+1-win)
        imax = min(len(L),i+1)
        Lf.append(sum(L[j] for j in range(imin, imax))/(imax-imin))
    plt.figure(4)
    plt.plot(d, L)
    plt.plot(d, Lf,lw=5,alpha=0.8)
    plt.show(block=False)

    return Lf


def ecdf(data, normalize=-999.0): # Incorporating Large-scale fading into the prediction
    '''
    generates CDF from empirical data 
    data - empirical data as a numpy vector
    normalize - percentil (%) to set zero (no normalization if outside <0,100>)
    returns - (x, y) vectors in tuple    
    '''
    x = np.sort(data)
    if 0 <= normalize <= 100.0: 
        x -= np.percentile(data, normalize)
    y = np.arange(1, len(x)+1)/float(len(x))
    return (x, y)


def one_slope_fit(d_cut, L_cut, Meas_name, plot, index):
    L_cut = np.array(L_cut)
    logd = np.log10(d_cut)
    n, L1, rv, pv, stderr = sst.linregress(logd, L_cut) # fit the 
    Lp = L1 + n * logd # prediction by the derived model

    err = Lp-L_cut
    stdev = np.std(err)
    
    if plot: 
         # show plot
        plt.figure(figsize=(9, 5))
        plt.semilogx(d_cut, L_cut, '+', lw=0.1, label='measurement')
        plt.semilogx(d_cut, Lp, lw=3, label=f'1SM: y intercept = {L1:.1f} dB, n = {0.1*n:.1f}')
        plt.xlabel('distance')
        plt.ylabel('path loss (dB)')
        plt.title('Path loss prediction vs. measurement')
        plt.legend()
        plt.grid()
        plt.show(block=False)
        if not (index == None):
            plt.savefig(f"figures/{Meas_name}/1_slope/Prediction_one_slope_{index}", dpi= 500)
        else:
            plt.savefig(f"figures/{Meas_name}/1_slope/Prediction_one_slope", dpi= 500)

        # calculate prediction error and print results
        err = Lp-L_cut
        stdev = np.std(err)
        print('Empirical One Slope Model parameters:')
        print('L1 =', L1.round(1), 'dB')
        print('n =', (0.1*n).round(1))
        print('standard deviation =', stdev.round(1), 'dB')
        print('mean error =', np.mean(np.abs(err)).round(1), 'dB')          
        
        
        # Large-scale Fading (Shadowing) Statistics
        # PDF
        plt.figure()
        plt.hist(err, bins=40, density=True, label='Empirical PDF')
        # comparison with Log-Normal distribution
        pdfx = np.linspace(-30,30)
        ppdf = sst.norm.pdf(pdfx, scale=stdev)
        plt.plot(pdfx, ppdf, label='Normal distribution fit')
        plt.legend()
        plt.grid()
        plt.show(block=False)
        if not (index == None):
            plt.savefig(f"figures/{Meas_name}/1_slope/PDF_{index}", dpi= 500)
        else:
            plt.savefig(f"figures/{Meas_name}/1_slope/PDF", dpi= 500)

        plt.figure(figsize=(9, 7))
        xcdf, ycdf = ecdf(err, 50)
        plt.semilogy(xcdf, 100*ycdf, label='Empirical CDF')
        plt.semilogy(xcdf,100*sst.norm.cdf(xcdf,scale=stdev), label='Normal CDF fit')
        plt.ylim(0.1,100)
        plt.xlabel('Normalized predicted RSS (dB)')
        plt.ylabel('Percentage of locations RSS < predicted RSS (%)')
        plt.legend()
        plt.grid()
        plt.show(block=False)
        if not (index == None):
            plt.savefig(f"figures/{Meas_name}/1_slope/CDF_{index}", dpi= 500)
        else:
            plt.savefig(f"figures/{Meas_name}/1_slope/CDF", dpi= 500)

        Pn = 20 # Normalizing the power???
        plt.figure(figsize=(9, 9))
        plt.subplot(211)
        plt.semilogx(d_cut, Pn-L_cut, '+')
        plt.plot(d_cut, Pn-Lp, lw =4, label='50 % location exceeding this' )
        plt.plot(d_cut, Pn-Lp-9.5, lw =4, label='90 % location exceeding this' )
        plt.plot(d_cut, Pn-Lp+9.5, lw =4, label='10 % location exceeding this' )
        plt.xlabel('Distance (m)')
        plt.ylabel('RSS (dBm)')
        plt.title('Predicted RSS (Distance in log and linear scale)')
        plt.legend()
        plt.grid()
        plt.subplot(212)
        plt.plot(d_cut, Pn-L_cut, '+')
        plt.plot(d_cut, Pn-Lp, lw =4, label='50 % location exceeding this' )
        plt.plot(d_cut, Pn-Lp-9.5, lw =4, label='90 % location exceeding this' )
        plt.plot(d_cut, Pn-Lp+9.5, lw =4, label='10 % location exceeding this' )
        plt.xlabel('Distance (m)')
        plt.ylabel('RSS (dBm)')
        plt.grid()
        if not (index == None):
            plt.savefig(f"figures/{Meas_name}/1_slope/RSS_{index}", dpi= 500)
        else:
            plt.savefig(f"figures/{Meas_name}/1_slope/RSS", dpi= 500)
        plt.show(block=True)

    # fit_matrix is the matrix with all the data that we got from the calculation
    return([Lp, L1, n, err, stdev]) # cislo rezu, odhadnute n, utlm v 1 m z prelozenia, predikcia utlmu, chyba predikcie, stdev


def cutout(d, d_start, d_end, L, cut_start, cut_end):
    # orezanie dlzok d a L
    if cut_start >= d_start and cut_end <= d_end:
        # sam_d = int(len(d)/(d_end-d_start))
        # sam_L = int(len(L)/(d_end-d_start))

        # s_in = (cut_start - 1)*sam_d # approximation of the index where the cutout STARTS
        # e_in = cut_end*sam_d - 1 # approximation of the index where the cutout ENDS

        s_in = 0
        e_in = 0
        
        # find the exact index
        # the int() approximation in lines 182 and 183 makes the calculated index smaller than it should have been
        # therefore the while cycles increase the index
        while d[s_in] < cut_start:
            s_in = s_in + 1
        e_in = s_in
        while d[e_in] < cut_end:
            e_in = e_in + 1
            if e_in == len(d):
                break
        
        d = d[s_in:e_in]
        L = L[s_in:e_in]
       
        return d, L
    else:
        print("ERR in the cutout function: Check your cut[] points")
        return d, L     


def slopes_fitting(d, d_start, d_end, L, cut, Meas_name): # Slope Model fitting

    if len(cut) > 0:
        d_matrix = []
        L_matrix = []
        L_model_matrix = [['Lp', 'L1', 'n', 'err', 'stdev']]        
            
        for c in range(1, len(cut)):
            d_cut, L_cut = cutout(d, d_start, d_end, L, cut[c-1], cut[c])
            d_matrix.append(d_cut)
            L_matrix.append(L_cut)
            
            print('-----------------------------------------------------------')
            print(f'Cutout interval number {c}: from {cut[c-1]} to {cut[c]} m :')
            print('[')
            L_model_matrix.append(one_slope_fit(d_cut, L_cut, Meas_name, True, c))
            print(']')


        return d_matrix, L_matrix, L_model_matrix
    
    else:
        return [], [], one_slope_fit(d,L)


if __name__ == "__main__":
    Main()
