%--------------------------------------------------------------------------
% This is the main function to run the SSC algorithm for the motion
% segmentation problem on the Hopkins 155 dataset.
%
% cd to the main folder containing the Hopkins 155 sequences
% add the path to the folder "SSC_motion_face" containing these m-files
%
% avgmissrate1: the n-th element contains the average clustering error for 
% sequences with n motions (using 2F-dimensional data)
% avgmissrate2: the n-th element contains the average clustering error for 
% sequences with n motions (using 4n-dimensional data)
% medmissrate1: the n-th element contains the median clustering error for 
% sequences with n motions (using 2F-dimensional data)
% medmissrate2: the n-th element contains the median clustering error for 
% sequences with n motions (using 4n-dimensional data)
%--------------------------------------------------------------------------
% Copyright @ Ehsan Elhamifar, 2012
%--------------------------------------------------------------------------

clc, clear all, close all

cd '/Users/ehsanelhamifar/Documents/MatlabCode/Hopkins155/';
addpath '/Users/ehsanelhamifar/Documents/MatlabCode/SSC_motion_face/';

alpha = 800;
    
maxNumGroup = 5;
for i = 1:maxNumGroup
    num(i) = 0;
end

d = dir;
for i = 1:length(d)
    if ( (d(i).isdir == 1) && ~strcmp(d(i).name,'.') && ~strcmp(d(i).name,'..') )
        filepath = d(i).name;
        eval(['cd ' filepath]);
        
        f = dir;
        foundValidData = false;
        for j = 1:length(f)
            if ( ~isempty(strfind(f(j).name,'_truth.mat')) )
                ind = j;
                foundValidData = true;
                break
            end
        end
        eval(['load ' f(ind).name]);
        cd ..
        
        if (foundValidData)
            n = max(s);
            N = size(x,2);
            F = size(x,3);
            D = 2*F;
            X = reshape(permute(x(1:2,:,:),[1 3 2]),D,N);
            
            r = 0; affine = true; outlier = false; rho = 0.7;
            [missrate1,C1] = SSC(X,r,affine,alpha,outlier,rho,s);
            
            r = 4*n; affine = true; outlier = false; rho = 0.7;
            [missrate2,C2] = SSC(X,r,affine,alpha,outlier,rho,s);

            num(n) = num(n) + 1;
            missrateTot1{n}(num(n)) = missrate1;
            missrateTot2{n}(num(n)) = missrate2;
            
            eval(['cd ' filepath]);
            save SSC_MS.mat missrate1 missrate2 C1 C2 alpha
            cd ..
        end   
    end
end

L = [2 3];
for i = 1:length(L)
    j = L(i);
    avgmissrate1(j) = mean(missrateTot1{j});
    medmissrate1(j) = median(missrateTot1{j});
    avgmissrate2(j) = mean(missrateTot2{j});
    medmissrate2(j) = median(missrateTot2{j});
end
save SSC_MS.mat missrateTot1 avgmissrate1 medmissrate1 missrateTot2 avgmissrate2 medmissrate2 alpha