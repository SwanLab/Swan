
clear;
close all;
clc;

etaVec = [0.5, 1.5];

for i = 1:length(etaVec)
    GrippingDensityCoupled(etaVec(i));
    GrippingLevelSetCoupled(etaVec(i));
end