function [ norma ] = normL2( unitM,fun )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
n = size(fun,2);
norma = 0;
for i = 1: n
    norma = fun(:,i)'*unitM*fun(:,i) + norma;
end
norma = sqrt(norma);
end

