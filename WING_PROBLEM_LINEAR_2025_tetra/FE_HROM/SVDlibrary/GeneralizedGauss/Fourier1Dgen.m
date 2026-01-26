function [f,df,DATAOUT] = Fourier1Dgen(xLIMloc,x,DATALOC)
% See Fourier1Dgen_aux.mlx
if nargin == 0
    load('tmp.mat')
end
PLOT_LOCAL =0;
%df = [] ;
DATAOUT = [] ;
N =DATALOC.PORDER ;  % Number of fourier functions
NFUN = 1+2*N ;
f = ones(length(x),NFUN) ;
df = zeros(length(x),NFUN) ;
acum = 1;
for n = 1:N
    acum = acum+1;
    f(:,acum) = cos(2*pi*x*n) ;
    df(:,acum) = -2*pi*n*sin(2*pi*x*n) ;
    acum = acum+1;
    f(:,acum) = sin(2*pi*x*n) ;
    df(:,acum) = 2*pi*n*cos(2*pi*x*n) ;
end

df = {df} ; 




if  PLOT_LOCAL == 1
    %  figure(45)
    
    
    xlabel('x')
    ylabel('f')
    hold on
    grid on
    for iplot = 1:size(f,2)
        plot(x,f(:,iplot)) ;
    end
end
