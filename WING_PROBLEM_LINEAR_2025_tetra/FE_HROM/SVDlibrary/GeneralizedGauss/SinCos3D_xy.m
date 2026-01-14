function   [A, dA,DATAOUT ]=    SinCos3D_xy(dummy,x,DATA)
%  f = sin(x*mu)*cos(y*mu)*z*mu
% ------------------------------------------------------------------------
% JAHO,  21-Oct-2021
% --------------------------____---------------------------------------------
if nargin == 0
    load('tmp2.mat')   
end
DATAOUT = [] ; 
nrows = size(x,1) ; 
ncols = length(DATA.MU) ; 
A = zeros(nrows,ncols) ; 
dA = cell(1,3) ; 
for i = 1:length(dA)
    dA{i} = zeros(nrows,ncols) ; 
end
for iparam = 1:length(DATA.MU)
    mu = DATA.MU(iparam) ; 
    A(:,iparam) =1+ sin(x(:,1)*mu).*cos(x(:,2)*mu).*(x(:,3)*mu) ; 
    idim = 1; 
    dA{idim}(:,iparam) =  mu*cos(x(:,1)*mu).*cos(x(:,2)*mu).*(x(:,3)*mu) ; 
      idim = 2; 
    dA{idim}(:,iparam) =  -mu*sin(x(:,1)*mu).*sin(x(:,2)*mu).*(x(:,3)*mu) ; 
     idim = 3; 
    dA{idim}(:,iparam) =  mu*sin(x(:,1)*mu).*cos(x(:,2)*mu)*mu ; 
    
end

% Constant function 
nA =  ones(size(A,1),1) ;
A = [A,nA] ; 
for idim = 1:3
    dA{idim} = [dA{idim},zeros(size(nA))] ; 
end


