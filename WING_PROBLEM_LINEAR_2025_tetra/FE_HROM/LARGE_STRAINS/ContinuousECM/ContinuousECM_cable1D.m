function [CECMoutput]= ContinuousECM_cable1D(BstRED_l,BasisTENSIONV,DATA,wFE,DATAoffline,DATA_ECM,...
    MESH,Nst)
% JAHO, 13-JAN-2021
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/ContinuousEmpiricalCubatureM/01_2D_beam_LARGE/README_OECM.pdf
% Adaptation to cables 1D :
%/home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/ContinuousEmpiricalCubatureM/Paper_hernandez2021ecm/09_ASSESSMENTpaper/FE1D_mooring/README_QuadMooring.mlx
% 31-OCT-2022


if nargin == 0
    load('tmp1.mat')
    DATAoffline.LoadFromMemory_ECMdata = 1; 
    DATAoffline.MAKE_VIDEO_POINTS = 1; %= DefaultField(DATAoffline,'MAKE_VIDEO_POINTS',0) ;
    DATAoffline.BasisIntegrandToPlot_INDEXES = [1:5,6:10:50] ; 
end
ECMdata = [] ;
DATA_ECM = DefaultField(DATA_ECM,'DIRECT_POLYNOMIAL_FITTING_INTEGRAND',1) ; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Matrix of reduced internal forces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A = BasisF_from_BasisStress_TENSION(BstRED_l,BasisTENSIONV,DATA)  ;

xNODES = MESH.COOR'; 
xFE = Nst*xNODES(:) ;
ndim = size(MESH.COOR,2) ;
xFE = xFE(1:ndim:end) ; 
MESH.COOR = MESH.COOR(:,1)  ; 

 
DATAMISC.BASE_FOLDER =   [cd,filesep,'CECMinfo',filesep];
if  ~exist( DATAMISC.BASE_FOLDER)
    mkdir( DATAMISC.BASE_FOLDER)
end
DATAMISC.ExactIntegral = A'*wFE ; 
DATA_ECM.TOL_SVD_A = DATAoffline.errorFINT ; 
MESH.ngausE = DATA.MESH.ngaus_RHS ; 

DATAoffline = DefaultField(DATAoffline,'LoadFromMemory_ECMdata',0) ; 
DATA_ECM = DefaultField(DATA_ECM,'NAMEws_CECMinfo',[cd,'/CECMinfo/','tmp_info.mat']) ; 


if DATAoffline.LoadFromMemory_ECMdata == 0
[CECMoutput,DATA_AUX]= ContinuousECMgen(A,xFE,wFE,DATAMISC,MESH,DATA_ECM) ;

save(DATA_ECM.NAMEws_CECMinfo,'CECMoutput','DATA_AUX') ; 

       

else
    load(DATA_ECM.NAMEws_CECMinfo,'CECMoutput','DATA_AUX') ;
    DATAoffline = DefaultField(DATAoffline,'MAKE_VIDEO_POINTS',0) ;
    DATA_AUX.AUXVAR = DefaultField(DATA_AUX.AUXVAR,'DATA_from_MAIN',[]) ;
    DATA_AUX.AUXVAR.DATA_from_MAIN.MAKE_VIDEO_POINTS = DATAoffline.MAKE_VIDEO_POINTS ; 
    DATAoffline = DefaultField(DATAoffline,'BasisIntegrandToPlot_INDEXES',[]) ; 
    DATA_AUX.AUXVAR.DATA_from_MAIN.BasisIntegrandToPlot_INDEXES = DATAoffline.BasisIntegrandToPlot_INDEXES ; 
  
end

VariousPLOTS_CECM(DATA_AUX.DATA_ECM,CECMoutput,DATA_AUX.MESH,DATA_AUX.AUXVAR,...
    DATA_AUX.VAR_SMOOTH_FE)  ;


        
[PHIk_y,~,~]=     EvaluateBasisFunctionALL(CECMoutput.xDECM,DATA_AUX.DATA_ECM,DATA_AUX.VAR_SMOOTH_FE,[]) ;
b=   DATA_AUX.DATA_ECM.ExactIntegral; % PHI'*W ;


bNEW = PHIk_y'*CECMoutput.wDECM ;
errorINT = norm(bNEW-b)/norm(b)*100;
disp(['Error DECM rule in integrating basis functions, including interpolation error (%) =',num2str(errorINT)]) ;
    
% 
    
    
[PHIk_y,~,~]=     EvaluateBasisFunctionALL(CECMoutput.xCECM,DATA_AUX.DATA_ECM,DATA_AUX.VAR_SMOOTH_FE,[]) ;
b=   DATA_AUX.DATA_ECM.ExactIntegral; % PHI'*W ;


bNEW = PHIk_y'*CECMoutput.wCECM ;
errorINT = norm(bNEW-b)/norm(b)*100;
disp(['Error CECM rule in integrating basis functions (%) =',num2str(errorINT)]) ;
    
% 
% %%%
% ECMdata.COOR = xNEW ;
% ECMdata.wRED = wNEW ;
% ECMdata.setPoints = []  ;
% ECMdata.setElements = ELEMENTS_xNEW ;
% ECMdata.Ninterpolation = Ninterpolation ;







