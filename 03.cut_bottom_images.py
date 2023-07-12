#!/usr/bin/python3
#
'''
Cut 31 px at the bottom of images downloaded from EcoTaxa.
Do it for both the learning set and the test set.
'''


import glob
import os
import lycon
import matplotlib.pyplot as plt
import numpy as np
import cv2

img_names = glob.glob('data/dataset/test/images/*/*.png')
n_img = len(img_names)

# Loop over images
for i, name in enumerate(img_names):

    # Read image
    img = lycon.load(name)[:,:,0]
    
    # Crop 31 px at the bottom
    h, w = img.shape
    img = img[0:h-31,:]
    
    # Get non blank parts and resize
    img_mask = img < 254
    h = np.max(np.where(img_mask)[0])
    w = np.max(np.where(img_mask)[1])
    img = img[0:h, 0:w]
    
    # Save image
    w = cv2.imwrite(name, img)
    
    if (i+1 % 1000 == 0):
        print(f'Done with {i+1} out of {n_img}')
    

