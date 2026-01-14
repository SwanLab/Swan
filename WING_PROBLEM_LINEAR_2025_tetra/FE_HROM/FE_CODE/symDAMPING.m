clc
clear all
% ---------------
%syms xiDmin alphaD betaD freqNmin freqNmax

%SOL = solve( xiDmin == 1/2*(alphaD/freqNmin + betaD*freqNmin), 1 == 1/2*(alphaD/freqNmax + betaD*freqNmax),alphaD,betaD)

freqNmin = 190 ; 
freqNmax = 5000 ; 
xiDmin = 0.01 ; 

 alphaD = (2*freqNmin*(- xiDmin*freqNmax^2 + freqNmin*freqNmax))/(freqNmin^2 - freqNmax^2) ; 
 betaD = -(2*(freqNmax - freqNmin*xiDmin))/(freqNmin^2 - freqNmax^2) ; 


% Least squares 
b = [xiDmin 1]' ; 
A = [0.5/freqNmin 0.5*freqNmin
     0.5/freqNmax  0.5*freqNmax]
 
SOL =  lsqnonneg(A,b)
xi = (A*SOL)

betaD = SOL(2); 
alphaD = SOL(1) ;

freqS = linspace(freqNmin,freqNmax,100) ; 
xi = 0.5*(alphaD./freqS + betaD*freqS); 

plot(freqS,xi)
