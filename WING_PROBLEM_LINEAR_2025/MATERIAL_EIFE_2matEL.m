function [PROPMAT,DATA ]= MATERIAL_EIFE_2matEL(DATA)
% -----------------------------------------------------------------------------------
% 3. Material data  (linear elasticity)
% -----------------------------------------------------------------------------------


DATA.TYPE_CONSTITUTIVE_MODEL_ALL = 'SMALL_STRAINS_ELASTIC' ;  
 
% EIF ELEMENT TYPE 1
imat =1 ; % Index material
NameEIFE_file = 'EIFE_LIBRARY/DEF_Q8_wing_1.mat' ;
load(NameEIFE_file,'EIFEoper') ;
PROPMAT(imat).EIFE_prop =  EIFEoper  ; %
PROPMAT(imat).CriterionChooseParentDomain =  'MAX_Q_11'  ;  % Criterion for choosing the connectivity of the element
PROPMAT(imat).PROPMAT = [] ;


PROPMATloc = [] ; 
jmat =1 ; % Index material
% Elasticity matrix
E = 70000  ; %  MPa, Young's modulus
nu = 0.3; % Poisson's coefficient
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G = E/2/(1+nu) ;  % Shear modulus
celasINV3D = [1/E  -nu/E -nu/E  0 0 0
    -nu/E  1/E  -nu/E  0 0 0
    -nu/E  -nu/E  1/E  0 0 0
    0      0    0  1/G   0 0
    0      0      0  0 1/G 0
    0      0      0  0  0  1/G] ;
ElasticityMatrix = inv(celasINV3D) ;
PROPMATloc(jmat).ElasticityMatrix =  ElasticityMatrix  ; %
PROPMATloc(jmat).Density = 2.7e-3 ; % MKg/m3


% Submaterial 2
jmat =2 ; % Index material
jREF = 1; 

PROPMATloc(jmat) = PROPMATloc(jREF) ; 



%%%%%%%

PROPMAT(imat).PROPMAT  = PROPMATloc ; 



%%%%%%%% EIF ELEMENT TYPE 2

imat = 2; 
IREF = 1; 

PROPMAT(imat).PROPMAT = PROPMAT(IREF).PROPMAT   ; 
PROPMAT(imat).EIFE_prop =  EIFEoper  ; %
PROPMAT(imat).CriterionChooseParentDomain =  'MAX_Q_33'  ;  % Criterion for choosing the connectivity of the element 
 