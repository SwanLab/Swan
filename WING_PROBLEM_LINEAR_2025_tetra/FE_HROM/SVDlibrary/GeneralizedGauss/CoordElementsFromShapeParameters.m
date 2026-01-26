function X = CoordElementsFromShapeParameters(DATAOUTlagrag,DATA)
% Parametrization of the coordinates of the 2D Lagrangian element of nodes
% [DATAOUTlagrag.xNODES,DATAOUTlagrag.yNODES], in terms of the "distorsion
% parameter" DATA.IRREGULAR_SHAPE_PARAMETRIZATION_ALPHA
% See CoordElementsFromShapeParameters_aux.mlx

if nargin== 0
    load('tmp.mat')
end

% 3) Parameterization in terms of alpha
% --------------------------------------


% 1) Coordinates of the domain (see numbering order in LagrangePolynomial2D_B_B_aux.mlx)
DATA = DefaultField(DATA,'IRREGULAR_SHAPE_PARAMETRIZATION_ALPHA',[]) ;
alpha = DATA.IRREGULAR_SHAPE_PARAMETRIZATION_ALPHA ;

if isempty(alpha)
    alpha = 0 ;
end
% % Numbering of Nshape ---
% [xNODES,yNODES] = meshgrid(DATAOUTlagrag.xNODES,DATAOUTlagrag.yNODES) ;
% xNODES = xNODES' ;
% yNODES = yNODES' ;
% % See picture in LagrangePolynomial2D_B_B_aux.mlx
% xNODES = xNODES(:) ;
% yNODES = yNODES(:) ;

Xbase = DATAOUTlagrag.COORnodesELEMENTS ;


if length(alpha)==1
    % Parameterization in terms of one single variable
    % ------------------------------------------------
    % Shear
    % --------------------------------------------------
    %  Points belonging to the
    Xloc= [-1, 1, +1-alpha, -1
        -1, -1, 1, 1  ] ;
    
    
    
elseif length(alpha)==4
    
    % alpha = [ANG1,ANG2,LENGTH1,LENGTH2]
    
    P1 = [-1,-1]' ;
    P2 = [+1,-1]' ;
    P4 = P1 + [cosd(alpha(1)),sind(alpha(1))]'*alpha(3) ;
    P3 = P2 + [cosd(alpha(2)),sind(alpha(2))]'*alpha(4) ;
    
    Xloc =[P1,P2,P3,P4] ;
    
    
    
else
    error('Option not implemented')
end

N =@(xi,eta) (0.25*[(1-xi).*(1-eta), (1+xi).*(1-eta), (1+xi).*(1+eta), (1-xi).*(1+eta) ]);
Xcoor = N(Xbase(1,:)',Xbase(2,:)')*Xloc(1,:)' ;
Ycoor = N(Xbase(1,:)',Xbase(2,:)')*Xloc(2,:)' ;
%     III = (xNODES==xLIM(1,1)) + (xNODES==xLIM(1,2)) +(yNODES==xLIM(2,1)) + (yNODES==xLIM(2,2)) ;
%     IND_BOUNDARIES =find(III>0) ;
%     XY_BND = [xNODES(IND_BOUNDARIES),yNODES(IND_BOUNDARIES)] ;
%     xy_BND = Ftransf*[XY_BND'] ;
PLOT_LOC =0;
if PLOT_LOC ==1
    figure(894)
    hold on
    xlabel('x')
    ylabel('y')
    plot([xNODES ],[yNODES ],'r*') ;
    axis equal
    plot([Xcoor  ],[Ycoor ],'bx') ;
    grid on
end
X = [Xcoor'; Ycoor'] ;
