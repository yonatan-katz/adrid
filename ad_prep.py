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
from scipy.fftpack import fft
import matplotlib.pyplot as plt

def make_filter(fc, b):        
   h = [
    0.000000000000000000,
    -0.001052835643627186,
    -0.004852904952195610,
    -0.002458804272575035,
    0.032350577876609214,
    0.118131300866432795,
    0.222462807119984562,
    0.270839718010742248,
    0.222462807119984646,
    0.118131300866432809,
    0.032350577876609235,
    -0.002458804272575040,
    -0.004852904952195614,
    -0.001052835643627185,
    0.000000000000000000,]    
   return h
    
# data = a numpy array containing the signal to be processed
# fs = a scalar which is the sampling frequency of the data   
def sfft(data, fs, fft_size=1024, overlap_fac=0.5): 
    hop_size = np.int32(np.floor(fft_size * (1-overlap_fac)))
    pad_end_size = fft_size          # the last segment can overlap the end of the data array by no more than one window size
    total_segments = np.int32(np.ceil(len(data) / np.float32(hop_size)))
    t_max = len(data) / np.float32(fs)
    T = 1 / fs
     
    window = np.hanning(fft_size)  # our half cosine window
    inner_pad = np.zeros(fft_size) # the zeros which will be used to double each segment size
     
    proc = np.concatenate((data, np.zeros(pad_end_size)))              # the data to process
    result = np.empty((total_segments, fft_size), dtype=np.float32)    # space to hold the result
    fig = 0
    for i in xrange(total_segments):                      # for each segment
        current_hop = hop_size * i                        # figure out the current segment offset
        segment = proc[current_hop:current_hop+fft_size]  # get the current segment
        windowed = segment * window                       # multiply by the half cosine function
        padded = np.append(windowed, inner_pad)           # add 0s to double the length of the data
        spectrum = np.fft.fft(padded)        
        
        xf = np.linspace(0.0, 1.0/(2.0*T), fft_size)            
        plt.figure(fig)
        plt.plot(xf,1.0/fft_size * np.abs(spectrum[0:fft_size])) 
        fig +=1
        
        result[i, :] = spectrum[:fft_size]               # append to the results array       
    
     
    result = 20*np.log10(result)          # scale to db
    result = np.clip(result, -40, 200)    # clip values
    
    return result


def get_sig_test(x, N=1024*4):    
    
    y1 = np.sin(10.0 * 2.0*np.pi*x[0:1024]) + 0.5*np.sin(150.0 * 2.0*np.pi*x[0:1024])    
    y2 = np.ones(N-1024)
    y2 = np.sin(20.0 * 2.0*np.pi*x[1024:]) + 0.5*np.sin(100.0 * 2.0*np.pi*x[1024:])    
    Y = np.append(y1,y2)    
    #y1 = np.sin(10.0 * 2.0*np.pi*x) +  \
    #    np.sin(150.0 * 2.0*np.pi*x)    
    return Y
    
def get_freq(sig):
    f = fft(sig)
    
    return f
    
def test(N=1024*4):
    T = 1.0 / 800.0    
    x = np.linspace(0.0, N*T, N)
    sig = get_sig_test(x,N)
    #f = make_filter(fc=0.1205,b=0.01)
    #sig = np.convolve(sig, f)    
    #fft = get_freq(sig)    
    #xf = np.linspace(0.0, 1.0/(2.0*T), N/2)    
    
    #plt.figure(0)
    #plt.plot(xf,2.0/N * np.abs(fft[0:N/2]))    
    
    spec = sfft(sig, fs=800.0)
    
    plt.figure(1)
    ##img = plt.imshow(spec, origin='lower', 
     #   cmap='jet', interpolation='nearest', 
     #   aspect='auto')
    
   # plt.show()
    
    return spec
    
    
    

def down_sample(signal):
    pass