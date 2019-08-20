function output = extract_TPM(input,choice,range)
%function output = extract_TPM(input,choice,range)

%This function calculates the Transition Probability matrices (TPM) of a
%given matrix along the 8 (N,S,W,E,NE,NW,SE,SW) directions for the 8
%graient matrices of the input matrix as explained in the paper
% Ravi, H., Subramanyam, A. V., Emmauel, S.,"Spatial Domain
% Quantization noise based Image Filtering Detection", in Proc. IEEE International 
% Conference on Image Processing (ICIP), Sep 2015. 
% Please cite the paper if you use this code.

%Inputs -
% input   - input M x N matrix for which the feature set has to be calculated. In the 
%	    paper it is 'noise' extracted using 'noise_finalone' code. 	
% choice  - 'sign'     - models the sign of the values in the matrix as states
%           'truncate' - quantizes and truncates the matrix between -range to +range and
%           then models the values itself as states.
% range   - is the -range to +range values to which the matrix will be
%           truncated if the choice is 'truncate'
%Outputs -
% output  - Feature set that concatenates TPM along wight directions.
%           Dimension length is ((2*range+1)^2)*8.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Code written by, @ Hareesh Ravi (Research Associate at IIITD)    %%%%  
%%%%  (haree.24@gmail.com)                                             %%%% 
%%%%  code can be used and modified for research purposes.             %%%% 
%%%%  Kindly let me know of mistakes by mail                           %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Difference array calculation
% Just two or four directions for certain applications might be optimal
dir{1} = input(2:end,:) - input(1:end-1,:); %Vertical bottom to top direction
dir{2} = input(1:end-1,:) - input(2:end,:); %Vertical top to bottom
dir{3} = input(:,2:end) - input(:,1:end-1); %Horizontal right to left 
dir{4} = input(:,1:end-1) - input(:,2:end); %Horizontal left to right
dir{5} = input(2:end,2:end) - input(1:end-1,1:end-1); %Major diagonal bottom to top
dir{6} = input(1:end-1,2:end) - input(2:end,1:end-1); %Minor diagonal top to bottom
dir{7} = input(2:end,1:end-1) - input(1:end-1,2:end); %Minor diagonal bottom to top
dir{8} = input(1:end-1,1:end-1) - input(2:end,2:end); %Major diagonal top to bottom

if strcmp(choice,'sign')
    val = 9;
elseif strcmp(choice,'truncate')
    val = (2*range+1)^2;
end

%TPM calculation
TPM = cell(1,8);
for i = 1 : 8
    D = dir{i};
    if strcmp(choice,'sign')
        D=sign(D);
        D=D+2;
    elseif strcmp(choice,'truncate')
        D=round(D);
        D(D>range)=range;D(D<-range)=-range;
        D=D+range+1;
    end
    %TPM calculates along the same direction that the gradient is
    %calculated (Can be modified if needed)
    TPM{i} = reshape(transprobestimate(D,choice,range,i),[1 val]);
end

output = [TPM{1},TPM{2},TPM{3},TPM{4},TPM{5},TPM{8},TPM{6},TPM{7}];
clear TPM;
end


function TPM = transprobestimate(input,choice,range,direction)
%function output = transprobestimate(input,choice,range,direction)

%Input -
% input     - input matrix
% choice    - 'sign' or 'truncate' as explained in 'extarct_TPM.m' function
% range     - range to be truncated 
% direction - values 1 to 8 for eight directions IN THE ORDER GIVEN IN
%             'extract_TPM.m' function.
%Output -
% TPM - TPM extracted from 'input', using 'choice' along 'direction'.

if strcmp(choice,'sign')
    t=zeros(3,3);
elseif strcmp(choice,'truncate')
    t=zeros(range*2+1,range*2+1);
end

%TPM calculation using inbuilt 'sparse' command. Better and optimized
%algorithms can easily be written without those commands. Hence though this
%code may not be optimized, it does the job. Here sparse command gives the
%occurences and converting it into a full matrix gives the co occurence
%matrix. Normalizing that, gives the probability matrix.

switch direction
    case 1
        %Vertical bottom to top direction
        t1=sparse(input(2:end,1:end),input(1:end-1,1:end),1);
    case 2
        %Vertical top to bottom direction
        t1=sparse(input(1:end-1,1:end),input(2:end,1:end),1);
    case 3
        %Horizontal right to left direction
        t1=sparse(input(1:end,2:end),input(1:end,1:end-1),1);
    case 4
        %Horizontal left to right direction
        t1=sparse(input(1:end,1:end-1),input(1:end,2:end),1);
    case 5
        %Major diagonal bottom to top direction
        t1=sparse(input(2:end,2:end),input(1:end-1,1:end-1),1);
    case 6
        %Minor diagonal top to bottom direction
        t1=sparse(input(1:end-1,2:end),input(2:end,1:end-1),1);
    case 7
        %Minor diagonal bottom to top direction
        t1=sparse(input(2:end,1:end-1),input(1:end-1,2:end),1);
    case 8
        %Major diagonal top to bottom direction
        t1=sparse(input(1:end-1,1:end-1),input(2:end,2:end),1);
    otherwise
        error('Wrong Argument direction');
end

t1=full(t1);
t(1:size(t1,1),1:size(t1,2))=t1;
for i=1:size(t,1)
    e(i,:) = t(i,:) ./ sum(t(i,:));
end
for i=1:size(e,1)
    for j=1:size(e,2)
        if (isnan(e(i,j)))
            e(i,j)=0;
        end
    end
end

TPM=e;
end