
clear;
close all;
clc;

etaVec = [3, 1/3];

for i = 1:length(etaVec)
    GrippingDensityCoupled(etaVec(i));
    GrippingLevelSetCoupled(etaVec(i));
end