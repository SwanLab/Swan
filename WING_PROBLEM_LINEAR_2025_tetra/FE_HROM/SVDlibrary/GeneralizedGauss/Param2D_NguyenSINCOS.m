function [f,df,DATAOUT] = Param2D_NguyenSINCOS(xLIMloc,x,DATALOC)
%   See Param2D_NguyenSINCOS_aux.mlx 
if nargin == 0
    load('tmp.mat')
  %  DATALOC.NumberOfSnapshots = [] ; 
end
%PLOT_LOCAL =0;
%df = [] ;
DATAOUT = [] ;
muDOMAIN =DATALOC.mu_ParametricDomain ;  % Parametric domain 
idim = 1; 
muGRID{idim} = linspace(muDOMAIN{idim}(1),muDOMAIN{idim}(2),DATALOC.NumberOfSnapshots(idim))' ;
idim = 2; 
muGRID{idim} = linspace(muDOMAIN{idim}(1),muDOMAIN{idim}(2),DATALOC.NumberOfSnapshots(idim))' ;

[xxMU,yyMU] = meshgrid(muGRID{1} ,muGRID{2} )  ; 

muSAMPLE = [xxMU(:),yyMU(:)] ; 


nsnapT = size(muSAMPLE,1) ; 

% include the function 1 
f = ones(size(x,1),nsnapT+1) ;
df = {zeros(size(f)),zeros(size(f))} ; 
 
 
for ifun = 1:nsnapT
    muLOC = muSAMPLE(ifun,:);
    idim = 1; 
    fx = sin(x(:,idim)*muLOC(idim)) ; 
    dfx = muLOC(idim)*cos(x(:,idim)*muLOC(idim)) ; 
    idim = 2; 
    fy = cos(x(:,idim)*muLOC(idim)) ; 
    dfy =  -muLOC(idim)*sin(x(:,idim)*muLOC(idim)) ; 
    
    f(:,ifun+1) = fx.*fy ; 
     
    df{1}(:,ifun+1) = dfx.*fy ; 
    df{2}(:,ifun+1) = fx.*dfy ;
  
end

%df = {df,df} ; 


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
