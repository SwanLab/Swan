%--------------------------------------------------------------------------
% This is the main function to run the SSC algorithm for the face 
% clustering problem on the Extended Yale B dataset.
% avgmissrate: the n-th element contains the average clustering error for 
% n subjects
% medmissrate: the n-th element contains the median clustering error for 
% n subjects
% nSet: the set of the different number of subjects to run the algorithm on
%--------------------------------------------------------------------------
% Copyright @ Ehsan Elhamifar, 2012
%--------------------------------------------------------------------------

clear all, close all

load YaleBCrop025.mat

alpha = 20;

nSet = [2 3 5 8 10];
for i = 1:length(nSet)
    n = nSet(i);
    idx = Ind{n};   
    for j = 1:size(idx,1)
        X = [];
        for p = 1:n
            X = [X Y(:,:,idx(j,p))];
        end
        [D,N] = size(X);
       
        r = 0; affine = false; outlier = true; rho = 1;
        [missrate,C] = SSC(X,r,affine,alpha,outlier,rho,s{n});
        missrateTot{n}(j) = missrate;
        
        save SSC_Faces.mat missrateTot alpha
    end
    avgmissrate(n) = mean(missrateTot{n});
    medmissrate(n) = median(missrateTot{n});
    
    save SSC_Faces.mat missrateTot avgmissrate medmissrate alpha
end

save SSC_Faces.mat missrateTot avgmissrate medmissrate alpha