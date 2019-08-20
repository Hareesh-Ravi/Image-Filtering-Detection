function [denoisedimage,noise,error]=noise_finalone(image,quality,weight,maxiter)
%function [denoisedimage,noise,error]=noise_finalone(image,quality,weight,maxiter)

% This algorithm removes compression noise from a given grayscale image
% whose quality factor is known. Based on the quality factor,a quantization
% noise model is created as given in the paper "DCT quantization noise in
% compressed images" by "Mark.A.Robertson and Robert.L.Stevenson". This
% code uses a modified HMRF prior model as discussed in our paper.
% Ravi, H., Subramanyam, A. V., Emmauel, S.,"Spatial Domain
% Quantization noise based Image Filtering Detection", in Proc. IEEE International 
% Conference on Image Processing (ICIP), Sep 2015.  
% Please cite the paper if you use this code.

%Input - 
% Image     - Noisy Image (mxn) - (0 to 255 range gray scale image)
% weight    - is the Huber function weight included as part of our application
%             (default is 1)
% maxiter   - number of iterations of gradient descent for every block
%             (default  - 10)
% quality   - QF of Jpeg compression used

%Output -
% denoisedimage - denoised Image (parameters have to be modified to get
%                 visual clarity. The default parameters are to get a good 
%                 estimate of compression noise and not visual clarity)
% noise - Compression noise extracted from compressed and denoised image.
% error - is the gradient descent error to check for convergence

%This code is to implement what has been explained in that paper for the
%specific problem. This code might not be optimized.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Code written by, @ Hareesh Ravi (Research Associate at IIITD)    %%%%  
%%%%  (haree.24@gmail.com)                                             %%%% 
%%%%  code can be used and modified for research purposes.             %%%% 
%%%%  Kindly let me know of mistakes by mail                           %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('weight','var')
    weight = 1;
end
if ~exist('maxiter','var')
    maxiter = 50;
end

%patchsize (DCT block size)
step = 8;

%Get the image
readcompressimage=double(image);

%initialize gamma as in paper
gamma = 1;

%Block size as givenby patchsize
blockX = step;blockY = step;

%Quantization table

Q50=[ 16 11 10 16 24 40 51 61;
     12 12 14 19 26 58 60 55;
     14 13 16 24 40 57 69 56;
     14 17 22 29 51 87 80 62; 
     18 22 37 56 68 109 103 77;
     24 35 55 64 81 104 113 92;
     49 64 78 87 103 121 120 101;
     72 92 95 98 112 100 103 99];

%Linear operation as done in JPEG compression
if quality > 50
    Q = round(Q50.*(ones(8)*((100-quality)/50)));
elseif quality < 50
    Q = round(Q50.*(ones(8)*(50/quality)));
elseif quality == 50
    Q = Q50;
end

    Q=double(Q);

% if step==16
%     Q = repmat(Q,2,2);
% end
% Create DCT Matrix
H = dct(eye(step));

% kronecker prod of DCT matrix
G =kron(H,H);

% getting noise covariance matrix
q=reshape(Q,[1 64]);

%determining variance as in the paper
varfreq = ((gamma.*q).^2)/12;
B = repmat(varfreq,[1 blockX*blockY/64]);
Key = double(diag(B));

%The noise covariance matrix
Kez = G'*Key*G;
kezinv=inv(Kez);

%initializing error matrix for grad descent. 
% k=1;error=zeros((size(readcompressimage,1)/step)*(size(readcompressimage,2)/step),maxiter);
%This can be uncommented and variable "error" canbe used inside the loop to
%get the error values to check for convergence

% Block by block denoising
for i=1:step:size(readcompressimage,1)
    for j=1:step:size(readcompressimage,2)
          
        y = readcompressimage(i:i+step-1,j:j+step-1);
%         zq=reshape(y,[step*step 1]);
        xblock = padarray(y,[1 1],'replicate');
        denoisedimage(i:i+step-1,j:j+step-1)= grad_finalone(y(:),xblock,kezinv,weight,maxiter);
%         k=k+1;
    end
end

%Extract noise from denoised image
noise = double(denoisedimage) - double(readcompressimage(1:size(denoisedimage,1),1:size(denoisedimage,2)));

noise = reshape(noise,[1 size(noise,1)*size(noise,2)]);
%truncating denoised image
denoisedimage(denoisedimage<0)=0;denoisedimage(denoisedimage>255)=255;


end

