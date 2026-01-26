function X = CoordElementsFromShapeParameters3D(DATAOUTlagrag,DATA)
% Parametrization of the coordinates of the 3D Lagrangian element of nodes
% [DATAOUTlagrag.xNODES,DATAOUTlagrag.yNODES], in terms of the "distorsion
% parameter" DATA.IRREGULAR_SHAPE_PARAMETRIZATION_ALPHA
% See CoordElementsFromShapeParameters3D_aux.mlx

if nargin== 0
    load('tmp.mat')
    DATA.IRREGULAR_SHAPE_PARAMETRIZATION_ALPHA=0.9 ; 
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
   
    Xloc= [-1, 1, +1, -1   -1 + alpha, 1-alpha, +1-alpha, -1+alpha 
          -1, -1, 1,  1     -1+alpha, -1+alpha, 1-alpha, 1-alpha 
          -1  -1  -1  -1   +1  +1  +1  +1] ;
    
% See CoordElementsFromShapeParameters3D_aux.mlx

    
 
    
elseif length(alpha)==2 
    % See CoordElementsFromShapeParameters3D_aux.mlx

     Xloc= [-1, 1, +1, -1   -1 , 1  , +1 , -1  
          -1, -1, 1,  1     -1 , -1  , 1-alpha(1) , 1-alpha(2)  
          -1  -1  -1  -1   +1  +1  +1  +1] ;
    
else
    error('Option not implemented')
end

N =@(xi,eta,zeta)  (0.125*[(1-xi).*(1-eta).*(1-zeta), (1+xi).*(1-eta).*(1-zeta), (1+xi).*(1+eta).*(1-zeta), (1-xi).*(1+eta).*(1-zeta),...
    (1-xi).*(1-eta).*(1+zeta), (1+xi).*(1-eta).*(1+zeta), (1+xi).*(1+eta).*(1+zeta), (1-xi).*(1+eta).*(1+zeta) ]);


Xcoor = N(Xbase(1,:)',Xbase(2,:)',Xbase(3,:)')*Xloc(1,:)' ;
Ycoor = N(Xbase(1,:)',Xbase(2,:)',Xbase(3,:)')*Xloc(2,:)' ;
Zcoor = N(Xbase(1,:)',Xbase(2,:)',Xbase(3,:)')*Xloc(3,:)' ;




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
    zlabel('z')
 %   plot3([xNODES ],[yNODES ],[zNODES],'r*') ;
    axis equal
    plot3([Xcoor  ],[Ycoor ],[Zcoor],'bx') ;
    grid on
end
X = [Xcoor'; Ycoor'; Zcoor'] ;
