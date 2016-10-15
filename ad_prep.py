#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  8 11:51:00 2016

@author: yonic
"""
'''Assuming get 44,1Khz mono signal
   and downsampling it to 11025Hz
'''
import numpy as np
import pandas as pd
from scipy.fftpack import fft
import matplotlib.pyplot as plt
import scipy.io.wavfile
from collections import deque

def make_filter():        
   #rate 44100
   #cutoff 11025
   #transition band 4000
   h = [
    0.001495793423956665,
    -0.000000000000000002,
    -0.002353465375814965,
    0.000000000000000003,
    0.004891460620271658,
    -0.000000000000000005,
    -0.009679017580867512,
    0.000000000000000008,
    0.017572376364284409,
    -0.000000000000000011,
    -0.030236514723564058,
    0.000000000000000015,
    0.051951159985932534,
    -0.000000000000000017,
    -0.098651926584576957,
    0.000000000000000019,
    0.315366534670899867,
    0.499287198398956633,
    0.315366534670899867,
    0.000000000000000019,
    -0.098651926584576957,
    -0.000000000000000017,
    0.051951159985932555,
    0.000000000000000015,
    -0.030236514723564062,
    -0.000000000000000011,
    0.017572376364284402,
    0.000000000000000008,
    -0.009679017580867519,
    -0.000000000000000005,
    0.004891460620271663,
    0.000000000000000003,
    -0.002353465375814965,
    -0.000000000000000002,
    0.001495793423956665,
    ]
   return h
    
# data = a numpy array containing the signal to be processed
# fs = a scalar which is the sampling frequency of the data   
# callback is a function that called to process moving window
# spectrum
def sfft(data, fs, callback, 
         is_plot = False,
         fft_size = 1024, 
         overlap_fac = 0.5): 
    hop_size = np.int32(np.floor(fft_size * (1-overlap_fac)))
    pad_end_size = fft_size          # the last segment can overlap the end of the data array by no more than one window size
    total_segments = np.int32(np.ceil(len(data) / np.float32(hop_size)))    
    T = 1.0 / fs
     
    window = np.hanning(fft_size)  # our half cosine window
    inner_pad = np.zeros(fft_size) # the zeros which will be used to double each segment size
     
    proc = np.concatenate((data, np.zeros(pad_end_size)))              # the data to process
    fig = 0
    for i in xrange(total_segments):                      # for each segment
        current_hop = hop_size * i                        # figure out the current segment offset
        segment = proc[current_hop:current_hop+fft_size]  # get the current segment
        windowed = segment * window                       # multiply by the half cosine function
        padded = np.append(windowed, inner_pad)           # add 0s to double the length of the data
        spectrum = np.fft.fft(padded)        
        
        if is_plot:        
            xf = np.linspace(0.0, 1.0/(2.0*T), fft_size)            
            plt.figure(fig)
            plt.plot(xf,1.0/fft_size * np.abs(spectrum[0:fft_size])) 
            fig +=1
            
        if callback is not None:
            xf = np.linspace(0.0, 1.0/(2.0*T), fft_size)
            callback(i, xf, 1.0/fft_size * np.abs(spectrum[0:fft_size]))      
        
     
    
    
def moving_fft_handler(offset, hz, power):
    print offset
    print hz.shape
    print power.shape
    pass


class PeakMaker:    
    
    def __init__(this, window=6):
        this.queue = deque()
        this.buffer = []
        this.window=window
        this.sum = 0.0
        this.is_init = False
        this.peaks_time = []
        this.peaks_hz = []

    def handler(this, offset, hz, power):
        if not this.is_init:                        
            this.queue.append(np.sum(power))
            
            if len(this.queue) == this.window:
                this.is_init = True
                this.sum = np.sum(this.queue)
        
        else:
            s = np.sum(power)
            this.sum -= this.queue.popleft()
            this.sum += s
            this.queue.append(s)
            mean = this.sum / (this.window * len(power))            
            
            for i in xrange(len(power)):
                if power[i] > mean * 3:
                    this.peaks_time.append(offset)
                    this.peaks_hz.append(hz[i])
                    
    def get_spectogramm(this):
        df = pd.DataFrame()
        df["hz"] = pd.Series(this.peaks_hz)
        df["time"]  =pd.Series(this.peaks_time)
        
        return df
        
        
    
    
def make_signature():
    rate,data = \
        scipy.io.wavfile.read(
          "/home/yonic/repos/adrid/data/milki_mono.wav")
    print rate, data.shape    
    data = data[0:44100*10]
    f = make_filter()
    data = np.convolve(data, f)
    data = data[0::4]
    rate = rate/4    
    peak_maker = PeakMaker(window=10)
    sfft(data, 
        fs=rate, 
        fft_size=1024, 
        callback=peak_maker.handler)
    
    return peak_maker.get_spectogramm()
    
    
def make_coordinates(sp,signal_id):
    coordintes = []
    i = 0
    for i in xrange(0,len(sp)-8):
        anchor_hz = sp.hz.iloc[i]
        anchor_time = sp.time.iloc[i]
        for j in range(i+3,i+3+5,1):
            p_hz = sp.hz.iloc[j]
            p_delta=np.abs(sp.time.iloc[j]-anchor_time)
            p = (anchor_hz,p_hz,p_delta)
            a = (anchor_time, signal_id)
            coordintes.append((p,a))
            
    return coordintes
            
      

            
        
        
    

       
    