
clear;
close all;
clc;

etaVec = [1/3, 3];

for i = 1:length(etaVec)
    GrippingDensityCoupled(etaVec(i));
end