#    RTL data capture using pyrlsdr.
#    Repeat the following 
#       Read the capture parameters from a param file 
#       Wait for old data file to disappear, 
#       Read samples from RTL, then write new file 
#       Repeat above for M different parameter sets (bands) if required
#
#    Sample usage "python rtl_cap.py 3"   to capture samples from 3 bands
#    Requires RTL2832 to be plugged in and pyrtlsdr to be installed
#
#    Data files will use 2 byte raw format, named 'datm.bin' m=1..M
#    Param files are text eg Fc1=436e6; Gain=5; Fs=2e6;    stored in
#    a text file called rtlparm.txt   where m is 1..M
#    Bill 2016 

from __future__ import division
from __future__ import print_function
from rtlsdr import *
import pylab as mpl
import time 
import numpy 
import os
import sys 

Nsamples = 1024*1024*4

# default values ... 
#Fc1 = 436e6
#Fs =  2e6 
#Gain = 4

def main():

    @limit_calls(2)
    def test_callback(samples, rtlsdr_obj):
        print('  in callback')
        print('  signal mean:', sum(samples)/len(samples))

    sdr = RtlSdr()

    #print('  sample rate: %0.6f MHz' % (sdr.rs/1e6))
    #print('  gain: %d dB' % sdr.gain)
    
    if len(sys.argv)>1:
        sin = str(sys.argv[1])
        Nband = int(sin)
    else:
        Nband = 2
    print('RTL capture with number of bands = ' + str(Nband))
    print('with ' + str(Nsamples) + ' samples')
    print('')
           
    while (1): 

       for m in range(1,Nband+1):
           
           dfname = 'dat'+'%d'%m +'.bin'
           while (os.path.isfile(dfname)==1):
               print(' wait for file ' + dfname + ' to disappear'); 
               time.sleep(4)


           pfname = 'rtlpar'+'%d'%m +'.txt'
           try:
               pars = ''
               with open(pfname,'r') as pfile:
                   pars = pfile.read()
           
               #try:
               #print(pars)
               exec pars in locals() 
           except:
               print('*** failed to setup RTL parameters for '+pfname)

           if len(pars)>1:
               sdr.fc = Fc
               sdr.rs = Fs
               sdr.gain = Gain
               print('Fc %0.3f MHz '  %(sdr.fc/1e6) + ' Fs %0.3f MHz' %(sdr.rs/1e6) + ' Gain %0.1f dB' %Gain)
               print('Reading samples...')
               samples = sdr.read_samples(Nsamples)
               samples = samples * 128
               print(' std %3.1f bits' % numpy.std(samples))

               array16 = numpy.array(samples.real, 'int16')
               print('write binary file ', dfname)
               newfile = open(dfname, 'wb')
               array16.tofile(newfile)
               array16 = numpy.array(samples.imag, 'int16')
               array16.tofile(newfile)
               newfile.close()
           else:
               print('parameter length = 0??')
               
           #time.sleep(1)


       

    

    print('Done\n')
    sdr.close()

if __name__ == '__main__':
    main()
