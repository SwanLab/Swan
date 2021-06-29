function [ rho ] = besselJroots(startFreq,P)
% rho = besselJroots(startFreq,P) returns a vector rho of size [P,1]
% constituted of roots of the Bessel function of first kind J_0. When
% startFreq is chosen as 0, rho contains the P first roots of J_0. When
% startFreq > 0, the function looks for roots around startFreq, sorted in
% ascending order according to their distance to startFreq. When there are
% more than P/2 roots between 0 and startFreq, the number p+ (resp. p-) of 
% roots in rho that are greater (resp. smaller) than startFreq satisfy
% - p+ + p- = P 
% - |p+ - p-| <= 1
% - p- >= p+
% If there are less than P/2 roots between roots between 0 and startFreq,
% then rho contains all roots between 0 and startFreq, and the rest of the
% roots are the first roots of J0 that are greater than startFreq. 

if nargin==0
    run('besselJrootsVal.m');
    % Validation script
else
    if startFreq==0
        f  = @(r) besselj(0,r);
        df = @(r) - besselj(1,r);
        % Initial guess using asymptotics of Bessel functions
        init_guess = (1:P)'*pi - pi/4;
        % Simple newton method 
        rho = newton(f,df,init_guess);
    else
        % We first look for closest root 'firstFreq' of J_0 near startFreq:
        f  = @(r) besselj(0,r);
        df = @(r) - besselj(1,r);
        nStart = round((startFreq + pi/4)/pi);
        startFreq = nStart*pi - pi/4;
        firstFreq = newton(f,df,startFreq);
        
        % Initial guess with proper repartition around firstFreq
        init_guess = (1:P)'*pi + max(firstFreq+pi/4 - (fix(P/2)*pi+3*pi/4),0)-pi/4;
        rho = newton(f,df,init_guess);
        
        % Re-order roots according to their distance to firstFreq
        [~,I] = sort(abs(rho-firstFreq));
        rho = rho(I);
        
    end
    rho = rho(:);
end
end

