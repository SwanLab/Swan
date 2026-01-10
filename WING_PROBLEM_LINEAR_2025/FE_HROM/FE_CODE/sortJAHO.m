function [Asort,Iforward,Iback] = sortJAHO(A,METHOD)
%  
%  A(Iforward) = Asort
%  Asort(Iback) = A; 
%
if nargin ==1
    METHOD = 'ascend' ; 
end

[Asort,Iforward] = sort(A,METHOD) ; 
Iback = zeros(size(Iforward)) ; 
Iback(Iforward) = 1:length(Iback) ;