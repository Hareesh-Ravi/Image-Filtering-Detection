# Image Filtering Detection
This repository has the code to process and extract features outlined in "Spatial domain quantization noise based image filtering detection". This paper addresses the problem of identifying forged or manipulated images that have been filtered using the added noise during compression process.  


## Reference
```
@INPROCEEDINGS{Ravi_2015_ICIP, 
author={H. {Ravi} and A. V. {Subramanyam} and S. {Emmanuel}}, 
booktitle={2015 IEEE International Conference on Image Processing (ICIP)}, 
title={Spatial domain quantization noise based image filtering detection}, 
year={2015}, 
pages={1180-1184}, 
doi={10.1109/ICIP.2015.7350986}, 
month={Sep.},}
```

## Instructions
------------------------------------------------------------------------------
Image = imread('ucid00001.tif'); 
imwrite(Image,'jff.jpg','Quality',70);
JpgImage = imread('jff.jpg');
[denoisedimage,noise]=noise_finalone(JpgImage,70,1,10);
Finalfeat = extract_TPM(noise,'truncate',15)
------------------------------------------------------------------------------
1. Main function is 'noise_finalone.m'
   This function extracts compression noise and hence the quality factor of the JPEG image
   which is denoised is assumed to be known. Each function has necessary comments and instructions
   on using them. One other function is called inside this function.

2. The values of threshols, maxiter, weight are all empirically determined and varies with application. 
   For eg, for general denoising purposes, weight 1 is better wheras to extarct effective compression noise
   weight 5 is better as mentioned in the following paper. Please cite the following paper if you use this.

3. "Final feat" gives you the final feature vector to be used for SVM training and testing.
   'extract_TPM.m' and the instructions to use the same are provided in the code as comments. The code extracts 
   only first order dependencies and it is not optimized. The code was written just using inbuilt commands.
   Two variants for this code are using 'sign' of the values as the states or 'the values' itself as states. 
   The range of the values required is user defined. It is advised to use a value such that 95% of the values 
   in the matrix are within the range [-value,+value].

'sign' choice was used in our previous paper titled
" Ravi, H., Subramanyam, A. V., Gupta, G. and Avinash Kumar, B., Compression noise Based Video Forgery Detection
  in Proc. IEEE International Conference on Image Processing (ICIP), Oct 2014, PP 5352-5356"

whereas 'truncate' choice was used in the paper
" Ravi, H., Subramanyam, A. V., Emmanuel, S., Spatial Domain Quantization noise based Image Filtering Detection, 
  in Proc. IEEE International Conference on Image Processing (ICIP), Sep 2015." 
 