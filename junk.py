#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 14 11:11:39 2016

@author: yonic
"""

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
import scipy.io.wavfile


def make_mono_and_record(
    post_fix=1,
    folder_songs="/home/yonic/repos/adrid/data/songs/",
    folder_records="/home/yonic/repos/adrid/data/records/") :
    
    fname=folder_songs+"ad_mono_"+str(post_fix)+".wav"
    rate,data = scipy.io.wavfile.read(fname)   
    mono = data[:,0]   
    
    for i in range(10):
        start=int(
            np.random.uniform(
              low=rate,
              high=len(mono)-rate*10,
              size=1))
        
        record = mono[start:start+rate*10]  
    
        print start      
    
        fname=folder_records+"ad_"+str(post_fix)+"_record_"+str(i)+".wav"
        scipy.io.wavfile.write(fname, rate, record)

    
    