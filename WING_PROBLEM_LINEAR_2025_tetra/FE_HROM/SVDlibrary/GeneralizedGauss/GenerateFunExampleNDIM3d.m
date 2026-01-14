function [Xf,W,mu,xMAT,MPOINTS,XY] = GenerateFunExampleNDIM3d(M,P,DATA,DATAFUN)
% Snapshot matrix Xf and integration weights for the
% TEsting function: page 12 of "A ‘best points’ interpolation method for efficient approximation
% of parametrized functions", by N.C. Nguyen et al, 2007
%---------------------------------------------------
% % G(x;mu) = (1-x)*cos(3*pi*mu(x+1))*e^((-1+x)*mu)
% Domain
if nargin == 0
    load('tmp1.mat')
end
xLIM = DATAFUN.xLIM;
Dmu = DATAFUN.muP ;
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
[xxGRID, yyGRID, zzGRID] = meshgrid(x{1},x{2},x{3}) ;
xx = xxGRID(:);
% Coordinate y
yy = yyGRID(:);
zz = zzGRID(:) ; 
% Therefore
xMAT = [xx yy zz] ;
MESH{1} = [xxGRID] ; 
MESH{2} = [yyGRID] ; 
 MESH{3} = [zzGRID] ;
% -----------------------------------------------





% Snapshot matrix


% Grid for the parametric space
mu = linspace(Dmu(1),Dmu(2),P) ;

if DATAFUN.REPEATMATRIX>0
    f = DATAFUN.REPEATMATRIX ;
else
    f = 1 ;
end

disp('Generating snapshot matrix ...')
if DATAFUN.TYPE(1) ==10
    Pdim = (P+1)^2 ;
else
    Pdim = P ; 
end 
disp(['Size = ',num2str(f*M(1)*M(2)*M(3)*Pdim*length(DATAFUN.TYPE)*8e-6),' Mbytes'])


M = length(W);

%if DATA.PARTITIONED_Xf == 0
Xf = [] ;
% XfDER_x = [] ;
% XfDER_y = [] ;
x = xMAT(:,1) ;
y = xMAT(:,2) ;
z = xMAT(:,3) ; 

% 
% DATA.DATAREMOVEPOINTS.HOLE.TYPE = 'POLYGONAL'; 
% DATA.DATAREMOVEPOINTS.HOLE.POINTS = [-0.6 -0.6
%                                      0.6 -0.6
%                                      0.6 0.6
%                                      -0.6 0.6]; 
DATA.DATAREMOVEPOINTS = DefaultField(DATA.DATAREMOVEPOINTS,'HOLE',[] ) ;



for idim = 1:length(DATAFUN.TYPE)
    % dbstop('33')
    
    %     XflocDERx = zeros(M,P) ;
    %     XflocDERy = zeros(M,P) ;
    
    if  DATAFUN.TYPE(idim) ==1
        Xfloc = SinCosExpFun(mu,x,y) ;
    elseif DATAFUN.TYPE(idim) ==10
        Xfloc = Poly2Dfun(P,x,y) ;
    elseif DATAFUN.TYPE(idim) ==11
      
        Xfloc =    LagrangePolynomial3D(xLIM,P,xMAT) ;
    elseif DATAFUN.TYPE(idim) ==12
        Xfloc =      LagrangePolynomial_ALLD_B_B(xLIM,P,xMAT,DATA) ;
    else
        error('Option not implemented')
    end
    Xf = [Xf Xfloc] ;
    %     XfDER_x =  [XfDER_x XflocDERx] ;
    %     XfDER_y =  [XfDER_y XflocDERy] ;
end

if ~isempty(DATA.DATAREMOVEPOINTS.HOLE)
    switch  DATA.DATAREMOVEPOINTS.HOLE.TYPE
        case 'POLYGONAL'
            POINTS_pol = DATA.DATAREMOVEPOINTS.HOLE.POINTS ;
            POINTS_pol = [POINTS_pol;POINTS_pol(1,:)] ;
            INPol = inpolygon(x,y,POINTS_pol(:,1),POINTS_pol(:,2)) ; 
            Xf(INPol,:) = 0 ; 
        otherwise 
            error('Option not implemented')
    end
end
 