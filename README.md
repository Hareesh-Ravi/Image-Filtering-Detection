# Image Filtering Detection
This repository has the code to process and extract features outlined in ICIP 2015 paper "Spatial domain quantization noise based image filtering detection" ( [PDF](https://ieeexplore.ieee.org/document/7350986) ) and ACM TOMM 2016 paper "Forensic Analysis of Linear and Nonlinear Image Filtering Using Quantization Noise" ( [PDF](https://dl.acm.org/citation.cfm?id=2857069) ). This paper addresses the problem of identifying forged or manipulated images that have been filtered using the added noise during compression process. If you use this code, please cite any of the following papers like given below.


## Reference
```
@article{Ravi:2016:FAL:2901366.2857069,
 author = {Ravi, Hareesh and Subramanyam, A. V. and Emmanuel, Sabu},
 title = {Forensic Analysis of Linear and Nonlinear Image Filtering Using Quantization Noise},
 journal = {ACM Trans. Multimedia Comput. Commun. Appl.},
 issue_date = {June 2016},
 volume = {12},
 number = {3},
 month = mar,
 year = {2016},
 issn = {1551-6857},
 pages = {39:1--39:23},
 articleno = {39},
 numpages = {23},
 url = {http://doi.acm.org/10.1145/2857069},
 doi = {10.1145/2857069},
 acmid = {2857069},
 publisher = {ACM}
} 
```
or
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
```
Image = imread('ucid00001.tif'); 
imwrite(Image,'jff.jpg','Quality',70);
JpgImage = imread('jff.jpg');
[denoisedimage,noise]=noise_finalone(JpgImage,70,1,10);
Finalfeat = extract_TPM(noise,'truncate',15)
```
1. Main function is ```noise_finalone.m```
   This function extracts compression noise and hence the quality factor of the JPEG image
   which is denoised is assumed to be known. 

2. The values of threshols, maxiter, weight are all empirically determined and varies with application. 
   For eg, for general denoising purposes, weight 1 is better wheras to extarct effective compression noise
   weight 5 is better as mentioned in the following paper.

3. "Final feat" returned by ```extract_TPM``` has the final feature vector to be used for SVM training and testing. 
   The code extracts only first order dependencies and it is not optimized. The code was written just using inbuilt commands.
   Two variants for this code are using 'sign' of the values as the states or 'the values' itself as states. 
   The range of the values required is user defined. It is advised to use a value such that 95% of the values 
   in the matrix are within the range [-value,+value].

```sign``` choice was used in our previous paper titled
" Ravi, H., Subramanyam, A. V., Gupta, G. and Avinash Kumar, B., Compression noise Based Video Forgery Detection
  in Proc. IEEE International Conference on Image Processing (ICIP), Oct 2014, PP 5352-5356"

```truncate``` choice was used in the paper
" Ravi, H., Subramanyam, A. V., Emmanuel, S., Spatial Domain Quantization noise based Image Filtering Detection, 
  in Proc. IEEE International Conference on Image Processing (ICIP), Sep 2015." 
 
