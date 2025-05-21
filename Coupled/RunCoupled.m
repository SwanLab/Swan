
clear;
close all;
clc;

etaVec = [0.2, 1.5];

for i = 1:length(etaVec)
    GrippingDensityCoupled(etaVec(i));
    GrippingLevelSetCoupled(etaVec(i));
end