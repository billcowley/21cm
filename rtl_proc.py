#    RTL data capture and processsing using pyrlsdr.
#    Repeat the following 
#       Read the capture parameters from a param file 
#       Wait for old PSD file to disappear, 
#       Read samples from RTL, find PSD, then write new file 
#       Repeat above for M different parameter sets (bands) if required
#
#    Sample usage "python rtl_proc.py"
#    Reads params from rtlpars.txt
#    Requires RTL2832 to be plugged in and pyrtlsdr to be installed
#
#    Bill 2018 

from __future__ import division
from __future__ import print_function
from rtlsdr import *
import pylab as mpl
from   scipy import signal
import time 
import numpy as np 
import os
import sys 

def main():

    @limit_calls(2)
    def test_callback(samples, rtlsdr_obj):
        print('  in callback')
        print('  signal mean:', sum(samples)/len(samples))

    sdr = RtlSdr()

    pfname = 'rtlpars.txt'
    try:
        pars = ''
        with open(pfname,'r') as pfile:
             pars = pfile.read()           
             exec pars in locals() 
    except:
        print('*** failed to setup RTL parameters for '+pfname)
        exit()    
    print('RTL capture with number of bands = ' + str(Nband))
    print('with ' + str(Nsam) + ' samples')
    print('')
    freqStored = 0
    
    while (1): 

       for m in range(1,Nband+1):
           
           dfname = 'psd'+'%d'%m +'.txt'
           while (os.path.isfile(dfname)==1):
               print(' wait for file ' + dfname + ' to disappear'); 
               time.sleep(4)
           try:        # reopen to allow changes on the fly 
               pars = ''
               with open(pfname,'r') as pfile:
                   pars = pfile.read()           
               exec pars in locals() 
           except:
               print('*** failed to setup RTL parameters for '+pfname)
               exit()

           if len(pars)>1:
               if (m==1): sdr.fc = Fc1
               if (m==2): sdr.fc = Fc2
               sdr.rs = Fs
               sdr.gain = Gain
               print('Fc %0.3f MHz '  %(sdr.fc/1e6) + ' Fs %0.3f MHz'
                     %(sdr.rs/1e6) + ' Gain %0.1f dB' %Gain)
               print('Reading samples...')
               samples = sdr.read_samples(Nsam)
               # returns floating point complex samples 
               print(' std %3.1f bits' % np.std(samples*128))

               # Welch PSD with default Hanning, no overlap
               f, Pxx = signal.welch(samples, fs=Fs, nperseg=Nw)
               print('PSD calculated')
               if (freqStored==0):
                   np.savetxt('freqs.txt', f, fmt='%12.1f')
                   freqStored=1
               Pxx= 10.0 * np.log10(Pxx)
               # write out the results in dBs    
               np.savetxt(dfname, Pxx, fmt = '%8.2f')
   
               
           else:
               print('bad params, skipping ??')
               
           #time.sleep(1)
           
    print('Done\n')
    sdr.close()

if __name__ == '__main__':
    main()
