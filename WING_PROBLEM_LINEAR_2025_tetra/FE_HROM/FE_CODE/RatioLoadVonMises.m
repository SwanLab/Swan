function ratioSTR = RatioLoadVonMises(nameMATERIAL,MATERIAL,stressVONMISES,DATA)

strenth = MATERIAL.(nameMATERIAL).strength ; 
stressCOMP = stressVONMISES.(nameMATERIAL) ; 
nfail = length(stressCOMP) ; 
nfail =  min(nfail,DATA.NUMBER_ELEMENTS_FAIL);% failed 
% Sort --> stressVONMISES.MATRIX
stressCOMP = sort(stressCOMP,'descend') ; 
% Maximum von mises stress for practical purposes 
stressMAX = stressCOMP(nfail); 
% Ratio between the Von Mises strength and the strength  
ratioSTR = stressMAX/strenth ; 