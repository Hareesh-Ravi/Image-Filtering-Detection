function [xopt,err] = grad_finalone(zq,xblock,kezinv,weight,maxiter)
%function [xopt,err] = grad_finalone(zq,xblock,z,kezinv,weight,maxiter)

%This function performs gradient descent for the gradient calculated from
%HMRF prior and quantization noise model as described in the paper
% Ravi, H., Subramanyam, A. V., Emmauel, S.,"Spatial Domain
% Quantization noise based Image Filtering Detection", in Proc. IEEE International 
% Conference on Image Processing (ICIP), Sep 2015.  
%Please cite the paper if you use this code.

%Input - 
% Zq - Noisy block vector (block^2 X 1) - (0 to 255 range gray scale block)
% xblock - padded noisy block of size (block+1 X block+1)
% weight - is the Huber function weight included as part of our application
% (default is 1)
% kezinv - inverse of the noise covariance matrix.
% maxiter - number of iterations of gradient descent for every block
% (default - 50)

%Output -
% xopt - Optimized denoised block after grad descent size (block X block)
% err - is the gradient descent error just to see if error is decreasing
%         with increasing number of iterations for every block

%This code is to implement what has been explained in that paper for the
%specific problem. This code might not be optimized.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Code written by, @ Hareesh Ravi (Research Associate at IIITD)    %%%%  
%%%%  (haree.24@gmail.com)                                             %%%% 
%%%%  code can be used and modified for research purposes.             %%%% 
%%%%  Kindly let me know of mistakes by mail                           %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


zq=double(zq);
z = zq;
% termination tolerance
tol = 1e-3;

% Initial step size for gradient descent
alpha = 0.05;

%Initial data vector
x = z;

%Initializing error vector
err=zeros(1,maxiter);
niter=1;

while (niter <= maxiter)
    
    
    % calculate gradient:
    g  = grad(xblock,z,zq,kezinv,weight);
    temp = z;

      % Gradient ddescent step
      z = z - alpha*g;
      
      %calculate error for tolerance measure
      err(1,niter) = mean(abs((abs(z) - abs(temp)))); 
   
      %upgrade alpha at every step
      alpha=alpha-(alpha*0.05);


      xnew = z;
      
      %error if value is Inf
      if ~isfinite(xnew)
            display(['Number of iterations: ' num2str(niter),'aplha is:'  num2str(alpha)])
            error('x is inf or NaN')
      end
    
      x = xnew;
      
      %Update xblock with new values to calculate gradient again
      [M,N] = size(xblock);
      xblockinter = reshape(xnew,M-2,N-2);

      xblock(2:end-1,2:end-1)=xblockinter;
      xblock(1,2:end-1)=xblock(2,2:end-1);
      xblock(end,2:end-1)=xblock(end,2:end-1);
      xblock(:,1)=xblock(:,2);
      xblock(:,end)=xblock(:,end-1);
      
      
      %come out of the loop if tol level reached
      if err(1,niter)<=tol
            break;
      end
      
      %Increment iteration counter
      niter = niter + 1;
end

% disp(niter-1);
%Final optimized/denoised block
[M,N] = size(xblock);
xopt = reshape(x,M-2,N-2);

end



function g  = grad(xblock,z,zq,kezinv,weight)
%function g  = grad(xblock,z,zq,kezinv,weight)

%This function implements the HMRF model & quantization noise model as 
%described in the same paper to a block and determines the gradient to be 
%used for the descent algorithm. 

% Input -
% xblock,z - block and vector of the block updated at every iteration
% zq - vector of noisy block
% kezinv - noise covariance matrix
% weight - introduced to govern the degree of denoising at the block
% boundaries

% Output -
% g - gradient calculated based on HMRF model and Quantization noise model.

%Lambda governs the degree of smoothing
lambda = 0.08;

%Threshold to define huber function as set in the paper
T = 10;
 
%Initializations
[M,N] = size(xblock);
xblock = double(xblock);
pho=zeros(3,3);%r=0;
phoo = zeros(5,5);
r=xblock(M-2,N-2);

%HMRF as in the paper
for l = 2:M-1
    for k = 2:N-1
        
        currentpx = xblock(l,k);
       
        for i = 1:3
            for j=1:3
                u = currentpx - xblock(l+i-2,k+j-2);   
                if(abs(u)<=T)
                    if(l-1==1 || k-1==1 || l-1==8 || k-1==8)
                       pho(i,j) = weight*(2*u);
                    else
                       pho(i,j) = 2*u;
                    end
                elseif((u)>T)
                    if(l-1==1 || k-1==1 || l-1==8 || k-1==8)
                        pho(i,j) = weight*2*T;
                    else
                        pho(i,j) = 2*T;
                    end
                elseif((u)<-T)
                    if(l-1==1 || k-1==1 || l-1==8 || k-1==8)
                        pho(i,j) = -weight*2*T;
                    else
                        pho(i,j) = -2*T;
                    end
                end

            end
        end
        
        phoo(2:4,2:4) = pho;
        for i = 2:size(phoo,1)-1
            for j=2:size(phoo,2)-1
                pho (i-1,j-1) = sum(sum(phoo(i,j) - phoo(i-1:i+1,j-1:j+1)));
            end
        end

        r(l-1,k-1) = sum(sum(pho));
    end
end

% 'r' term from HMRF
rr=reshape(r,[size(r,1)*size(r,2) 1]);

diffZ = double(z - zq);

% 's' term from quantization noise model
s = kezinv*diffZ; 

%gradient as computed in the paper
g = lambda.*rr + s;
end