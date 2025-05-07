
clear;
close all;
clc;

epsOverhVec = [2,4,6,8];
tarVec      = [0.01, 0.034, 0.067, 0.1];

for i = 1:length(epsOverhVec)
    for j = 1:length(tarVec)
        GrippingDensityLocalCirlceDomains(epsOverhVec(i),tarVec(j));
    end
end