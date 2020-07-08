function [ rho ] = besselYroots(startFreq,P)
% Same as besselJroots but replacing J_0 by Y_0 (Bessel function of second
% kind). See besselJroots for description. 

if nargin==0
    run('besselYrootsVal.m');
else
    if startFreq==0
        f  = @(r) bessely(0,r);
        df = @(r) - bessely(1,r);
        % Initial guess using asymptotics of Bessel functions
        init_guess = pi/4 + (0:P-1)'*pi;
        % Simple newton method
        rho = newton(f,df,init_guess);
    else
        % We first look for closest root 'firstFreq' of Y_0 near startFreq:
        f  = @(r) bessely(0,r);
        df = @(r) - bessely(1,r);
        nStart = round((startFreq - pi/4)/pi);
        startFreq = nStart*pi + pi/4;
        firstFreq = newton(f,df,startFreq);
        
        % Initial guess with proper repartition around firstFreq
        init_guess = (0:P-1)'*pi + max(firstFreq-pi/4 - (fix(P/2)*pi-pi/4),0)+pi/4;
        rho = newton(f,df,init_guess);
        
        % Re-order roots according to their distance to firstFreq
        [~,I] = sort(abs(rho-firstFreq));
        rho = rho(I);
        
    end
end
end

