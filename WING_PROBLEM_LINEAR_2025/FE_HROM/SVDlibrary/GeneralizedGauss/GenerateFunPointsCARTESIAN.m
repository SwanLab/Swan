function [Xf,W,xMAT,MPOINTS,XY] = GenerateFunPointsCARTESIAN(M,DATA)
% Snapshot matrix Xf and integration weights  
%---------------------------------------------------
if nargin == 0
    load('tmp2.mat')
end
xLIM = DATA.xLIM;


DATA = DefaultField(DATA,'INITIAL_DISCRETIZATION_TENSORPRODUCT',0);


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
    %MESH = cell(ndim,1) ; 
    %MESH{1} = [xxGRID] ;
    %MESH{2} = [yyGRID] ;
    %MESH{3} = [zzGRID] ;
    
elseif ndim==2 
      [xxGRID, yyGRID] = meshgrid(x{1},x{2}) ;
    xx = xxGRID(:);
    % Coordinate y
    yy = yyGRID(:);
    % Therefore
    xMAT = [xx yy  ] ;
    %MESH = cell(ndim,1) ; 
    %MESH{1} = [xxGRID] ;
    %MESH{2} = [yyGRID] ;
    
elseif ndim==1 
    xMAT =x{1} ; 
  
else
    error('Option Not implemented')
end
% -----------------------------------------------

M = length(W);

% Generation of snapshots 
% ---------------------------
% DATA.DATAFUN.TYPE =  'LagrangePolynomial3Dgen' ; % LagrangePolynomial_ALLD_B_Bgen % Function to be integrated 
DATA.EvaluateGradientsHere = 0; 
Xf = feval(DATA.TYPEFUN,xLIM,xMAT,DATA) ; 
% xLIM = Limits cartesian domain 
% xMAT = Points at which the functions are to be evaluated
 % DATA = Other inputs data relevant for the function "DATA.TYPEFUN"

 