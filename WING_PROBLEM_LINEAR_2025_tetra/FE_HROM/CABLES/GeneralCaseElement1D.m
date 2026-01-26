function MESH1D = GeneralCaseElement1D(DATAINPUT ) ;
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/ContinuousEmpiricalCubatureM/Paper_hernandez2021ecm/09_ASSESSMENTpaper/FE1D_mooring/README_QuadMooring.mlx
% Construction of meshes using elements of arbitrary order of interpolation
% JAHO, 31-Oct-2022
%----------------------------
if nargin == 0
    load('tmp1.mat')
end


% DATAINPUT.DATAcommon  = DefaultField(DATAINPUT.DATAcommon,'TypeElement_order_polynomial_1D',1)  ;
% DATAINPUT.DATAcommon  = DefaultField(DATAINPUT.DATAcommon,'Number_of_GaussPoints_per_element1D',[])  ;
%

nnodeE = DATAINPUT.TypeElement_order_polynomial_1D + 1 ;
ngausE = DATAINPUT.Number_of_GaussPoints_per_element1D ;

MESH1D.L =  DATAINPUT.LENGTH_CABLE; % m, initial length
NELEMS = DATAINPUT.NELEMS ;%
% 1 element --- n nodes
% 2 elements ---- 2*n-1
% 3 elements ---- 3*n-2
% m elements ---- m*n - (m-1) = m(n-1) + 1
NNODES = NELEMS*(nnodeE-1)+ 1;
COOR = zeros(NNODES,3) ;
COOR(:,1)  = linspace(0,MESH1D.L,NNODES);
MESH1D.COOR = COOR;

% CONNECTIVITIES
CN = zeros(NELEMS,nnodeE) ;
for inodeE = 1:nnodeE
    CN(:,inodeE)  = inodeE:(nnodeE-1):(NNODES-nnodeE+inodeE) ;
end
MESH1D.CN = CN;
MESH1D.TypeElement = 'Linear' ;

% Change the order for quadratic elements
if nnodeE == 3
    CNold = CN ;
    MESH1D.CN(:,2)  = CNold(:,3) ;
    MESH1D.CN(:,3)  = CNold(:,2) ;
end



MESH1D.POINTS_BOUNDARY = [1,NNODES] ;
MESH1D.MaterialType = ones(NELEMS,1) ;


if ~isempty(DATAINPUT.Number_of_GaussPoints_per_element1D)
    %  ngausLOC = [6,6];
    %  [~, ~, xGAUSS, weights]  =TensorProd2Ddiscr(ngausLOC) ;
    [posgp, weights] = GaussQuad(DATAINPUT.Number_of_GaussPoints_per_element1D,-1, +1) ;
    %     [xx,yy]  = meshgrid(xGAUSS{1},xGAUSS{2}) ;
    %     xx = xx(:) ; yy = yy(:);
    %     posgp= [xx'; yy'] ;
    MESH1D.posgp_given  = posgp(:)' ;
    MESH1D.weights_given  = weights ;
else
    %     MESH1D.posgp_given = [] ;
    %     MESH1D.weights_given = []  ;
    error('You should provide the number of Gauss points per element ')
end
