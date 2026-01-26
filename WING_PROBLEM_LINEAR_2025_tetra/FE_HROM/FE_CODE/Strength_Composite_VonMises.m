function  DATAOUT =   Strength_Composite_VonMises(DATAOUT,MATERIAL,DATA) 
 % ---------------------------------
 % Computing  ULTIMATE STRENGHT AND STRAINS (stress space and space
 % strains) . Von Mises failure criteria for both fiber and matrix
 % ---------------------------------   
 %dbstop('7')
 if nargin ==0
     load('tmp.mat')
 end
 
 % Principal stresses 
[ PRINCIPAL_STRESSES stressVONMISES]=  Principal_Stresses_Comp(DATAOUT,MATERIAL) ; 

% Number of elements that have to fail to consider that the RVE has also


ratioSTR =[] ; 

% MATRIX 
% .......
nameMATERIAL = 'MATRIX' ; 

ratioSTR_matrix = RatioLoadVonMises(nameMATERIAL,MATERIAL,stressVONMISES,DATA) ; 

nameMATERIAL = 'FIBER' ; 

ratioSTR_fiber = RatioLoadVonMises(nameMATERIAL,MATERIAL,stressVONMISES,DATA) ; 

% The   load factor to be employed is the maximum of ratioSTR_matrix and
% ratioSTR_fiber 
load_factor = max(ratioSTR_matrix,ratioSTR_fiber) ; 

% Finally, hence, the valu of the MACRO-stress vector when the RVE starts to fail
% is 
stressAVG_FAIL = DATAOUT.stressAVG/load_factor ; 
% Likewise, for strains 
strainAVG_FAIL = DATA.strainMACRO/load_factor ;  

DATAOUT.stressAVG_FAIL = stressAVG_FAIL ; 
DATAOUT.strainAVG_FAIL = strainAVG_FAIL ; 