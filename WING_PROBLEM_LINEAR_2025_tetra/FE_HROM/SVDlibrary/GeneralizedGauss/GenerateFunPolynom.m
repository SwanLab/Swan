function [Xf,W,xMAT,MPOINTS,XY] = GenerateFunPolynom(M,P,DATA,DATAFUN)
% Snapshot matrix Xf and integration weights for polynomials and 
% other used-defined functions 
%---------------------------------------------------
if nargin == 0
    load('tmp2.mat')
end
xLIM = DATAFUN.xLIM;

DATA = DefaultField(DATA,'PARTITIONED_Xf',0);
DATA = DefaultField(DATA,'PLOT_FUNCTION_SNAP',0);
DATAFUN = DefaultField(DATAFUN,'REPEATMATRIX',0);

% Number of integration points (spatial grid)
% Location of integration points (and associated weights)
DATA = DefaultField(DATA,'INITIAL_DISCRETIZATION_TENSORPRODUCT',0);
%dbstop('21')

% ----------------------------------------------------------
[MPOINTS, x, W] = InitalDiscretization(DATA,M,xLIM) ; 
% ndim = 2;
%dbstop('38')
XY = x; 
ndim = length(x) ; 
if ndim ==3
    [xxGRID, yyGRID, zzGRID] = meshgrid(x{1},x{2},x{3}) ;
    xx = xxGRID(:);
    % Coordinate y
    yy = yyGRID(:);
    zz = zzGRID(:) ;
    % Therefore
    xMAT = [xx yy zz] ;
    MESH = cell(ndim,1) ; 
    MESH{1} = [xxGRID] ;
    MESH{2} = [yyGRID] ;
    MESH{3} = [zzGRID] ;
    
else
    option('Not implemented')
end
% -----------------------------------------------

M = length(W);

% Generation of snapshots 
% ---------------------------
% DATA.DATAFUN.TYPE =  'LagrangePolynomial3D' ; % LagrangePolynomial_ALLD_B_B % Function to be integrated 

Xf = feval(DATA.DATAFUN.TYPE,xLIM,P,xMAT,DATA) ; 
% xLIM = Limits cartesian domain 
% P = Order of polynomial (or any other parameter)
% xMAT = Points at which the functions are to be evaluated 

 