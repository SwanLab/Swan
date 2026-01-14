clc
clear all
m = 1000 ;
n = 1000 ;
r = 100 ;
S1 = 10 ; 
mu = min(m,n)*eps ;
% We set lambda so that the last SV is equal to mu
lambda = -log(mu)/(r-1) ;
nmodes = 1:r ;

S =   S1*exp(lambda*(1-nmodes))';


SingVsq =  (S.*S) ;
SingVsq = sort(SingVsq);  % s_r, s_{r-1} ... s_1
normEf2 = sqrt(cumsum(SingVsq)) ; % s_r , s_r + s_{r-1} ... s_r +s_{r-1}+ s_{r-2} ... + s_1

figure(1)
hold on 
lS = log10(normEf2)

plot(lS)