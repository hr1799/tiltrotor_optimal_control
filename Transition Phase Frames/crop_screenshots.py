import cv2
import numpy as np

import os

files = [f for f in os.listdir('.') if os.path.isfile(f)]
count = 1
for f in files:
    if "Screenshot" in f:
        print(f)
        img = cv2.imread(f)

        marg_l = 0
        marg_r = 0
        marg_t = 300
        marg_b = 200
        
        # Cropping an image
        cropped_image = img[marg_t:1080-marg_b, marg_l:1920-marg_r, :]
        
        # Save the cropped image
        cv2.imwrite("frame"+str(count)+".png", cropped_image)
        count += 1