function Nshape = InverseMapping3Dinterpolation(MESHcoar,DATA,COORbnd)

if nargin == 0
    load('tmp1.mat')
end
% Quadratic elements require inverse mapping...
% OPTION INTRODUCED APRIL 3rd 2023
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/101_MULTI2D_2023/06_3D_27_HEXAHEDRA/04_FinalExam.mlx

TypeElement = MESHcoar.TypeElement;
DATA.InterpolationMethod = 'FEnodes'  ;
DATA = DefaultField(DATA,'ORDER_INVERSE_ELEMENT_inverse_mapping',2) ;
DATA = DefaultField(DATA,'NPOINTS_ONE_DIRECTION_inverse_mapping',10) ;


COORnodes = MESHcoar.COOR(MESHcoar.CN,:);

coorREF  = sum(COORnodes,1)/size(COORnodes,1) ;  %
%
COORnodes_REF = bsxfun(@minus,COORnodes',coorREF')';
COORbndREF = bsxfun(@minus,COORbnd',coorREF')';

LelemSCALING = max(abs(COORnodes_REF)) ;
%
% Transformed coordinates (reltive to coorREF), and scaled by   LelemSCALING
for idim = 1:length(LelemSCALING)
    COORnodes(:,idim) = COORnodes_REF(:,idim)/LelemSCALING(idim);
    COORbndREF(:,idim) = COORbndREF(:,idim)/LelemSCALING(idim);
end
%

[Nshape,~,~] = ShapeFunDer_inversemapping(COORnodes',COORbndREF,TypeElement,DATA) ;