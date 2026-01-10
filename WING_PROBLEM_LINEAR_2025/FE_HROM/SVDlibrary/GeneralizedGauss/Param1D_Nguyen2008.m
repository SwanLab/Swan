function [f,df,DATAOUT] = Param1D_Nguyen2008(xLIMloc,x,DATALOC)
% See Param1D_Nguyen2008_aux.mlx
if nargin == 0
    load('tmp.mat')
    DATALOC.NumberOfSnapshots = 4 ; 
end
PLOT_LOCAL =0;
%df = [] ;
DATAOUT = [] ;
muDOMAIN =DATALOC.mu_ParametricDomain ;  % Parametric domain 
nsnap = DATALOC.NumberOfSnapshots ; 
muSAMPLE = linspace(muDOMAIN(1),muDOMAIN(2),nsnap) ; 
% include the function 1 
f = ones(length(x),nsnap+1) ;
df = [];
 
for ifun = 1:nsnap
     muLOC = muSAMPLE(ifun); 
    f(:,ifun+1) = (1-x).*cos(3*pi*muLOC*(x+1)).*exp(-(1+x)*muLOC) ; 
  
end

df = {df} ; 




if  PLOT_LOCAL == 1
      figure(45)
    
    
    xlabel('x')
    ylabel('f')
    hold on
    grid on
    for iplot = 1:size(f,2)
        plot(x,f(:,iplot)) ;
    end
end
