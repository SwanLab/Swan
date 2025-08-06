
clear;
close all;
clc;

epsOverhVec = [4,6];
tarVec      = [0.034, 0.067];

for i = 1:length(epsOverhVec)
    for j = 1:length(tarVec)
        GrippingDensityLocalCirlceDomains(epsOverhVec(i),tarVec(j));
    end
end