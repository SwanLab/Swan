function [f,df,DATAOUT] = Param2D_Nguyen2008_poly_y(xLIMloc,x,DATALOC)
% SAme function as in Param1D_Nguyen2008 (see Param1D_Nguyen2008_aux.mlx),
% but along both x and y directions 
if nargin == 0
    load('tmp.mat')
  %  DATALOC.NumberOfSnapshots = [] ; 
end
%PLOT_LOCAL =0;
%df = [] ;
DATAOUT = [] ;
muDOMAIN =DATALOC.mu_ParametricDomain ;  % Parametric domain 
nsnap = DATALOC.NumberOfSnapshots ; 
muSAMPLE = linspace(muDOMAIN(1),muDOMAIN(2),nsnap)' ;

nsnapT = size(muSAMPLE,1) ; 

% include the function 1 
f = ones(size(x,1),nsnapT+1) ;
 
df = [];
 
for ifun = 1:nsnapT
    muLOC = muSAMPLE(ifun);
    idim = 1; 
    fx = (1-x(:,idim)).*cos(3*pi*muLOC*(x(:,idim)+1)).*exp(-(1+x(:,idim))*muLOC) ; 
%     idim = 2; 
%     fy = (1-x(:,idim)).*cos(3*pi*muLOC*(x(:,idim)+1)).*exp(-(1+x(:,idim))*muLOC) ;
    fy =  x(:,2); 
    
    f(:,ifun+1) = fx.*fy ; 
  
end

df = {df,df} ; 


% 
% 
% if  PLOT_LOCAL == 1
%       figure(45)
%     
%     
%     xlabel('x')
%     ylabel('f')
%     hold on
%     grid on
%     for iplot = 1:size(f,2)
%         plot(x,f(:,iplot)) ;
%     end
% end
