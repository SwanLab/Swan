function [PROPMAT, typePROBLEM, densGLO,DATA,CelasGLO] = ...
    AssignMatProp_J2(ndim,MATERIAL,nelem,MaterialType,DATA,TypeElement,nnodeE)

%dbstop('4')
if nargin == 0
    load('tmp.mat')
end

% Number of Gauss Points 
% ---------------------- 
[weig,posgp,shapef,dershapef] = ComputeElementShapeFun(TypeElement,nnodeE,'K') ;
ngausE = size(posgp,2) ; 

typePROBLEM = 'pstress';  %'pstress'/'pstrain'/'3D';  Plane stress/ plane strain problem
DATA = DefaultField(DATA,'StrainStressWith4Components',1) ;
if ndim==2 && DATA.StrainStressWith4Components == 0
    nstrain = 3;
elseif ndim==2 && DATA.StrainStressWith4Components == 1
    nstrain = 4 ;
    typePROBLEM = 'pstrain';
else
    nstrain = 6 ;
    typePROBLEM ='3D' ;
end

[DATA,PROPMAT,densGLO,CelasGLO] = J2initialVARIABLES(nstrain,DATA,nelem,MaterialType,MATERIAL,ngausE) ; 



%%%%%

%%%

