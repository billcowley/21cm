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
#
#   add option for discarding samples ... "python rtl_proc.py -d"
     
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
    sdrFail=0
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
    
    discardMode = 0
    if len(sys.argv)>1:
        if str(sys.argv[1])=="-d":
            discardMode = 1;
            print(" running in discard mode ")

    while (1): 

       for m in range(1,Nband+1):
           
           dfname = 'psd'+'%d'%m +'.txt'
           if discardMode==0: 
              while (os.path.isfile(dfname)==1):
                 #print(' wait for file ' + dfname + ' to disappear');
                 sys.stdout.write('w');   sys.stdout.flush()
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
               print('collect data from chan'+str(m))
#               print(sdr)
               try: 
                  if (m==1): sdr.fc = Fc1
                  if (m==2): sdr.fc = Fc2
                  sdr.rs = Fs
                  sdr.gain = Gain
                  print('Fc %0.3f MHz '  %(sdr.fc/1e6) + ' Fs %0.3f MHz'
                     %(sdr.rs/1e6) + ' Gain %0.1f dB' %Gain)
                  print('Reading samples...')
                  samples = sdr.read_samples(Nsam)
               except:
                  print("sdr failed")
                  sdrFail=1
                  time.sleep(1)
                  sdr.close()
                  time.sleep(1)
                  break
                        
               # demean has little effect          samples = samples - np.mean(samples)
               # returns floating point complex samples 
               print(' std %3.1f bits' % np.std(samples*128))

               # Welch PSD with default Hanning, no overlap
               f, Pxx = signal.welch(samples, fs=Fs, detrend = 'constant', nperseg=Nw)
               Pxx[0] = (Pxx[3]+Pxx[2]+Pxx[-2]+Pxx[-3])/4.0   # replace DC value
               Pxx[1] = Pxx[0];
               Pxx[-1]= Pxx[0];  

               print('PSD calculated')
               if (freqStored==0):
                   np.savetxt('freqs.txt', f, fmt='%12.1f')
                   freqStored=1
               Pxx= 10.0 * np.log10(Pxx)
               # write out the results in dBs    
               np.savetxt(dfname, Pxx, fmt = '%8.2f')
   
               
           else:
               print('bad params, skipping ??')
               
           time.sleep(0.1)
       if sdrFail==1:
           time.sleep(1)
           sdr = RtlSdr()
           sdrFail = 0 # hopefully reopen RTL dongle 
           
    print('Done\n')
    sdr.close()

if __name__ == '__main__':
    main()
