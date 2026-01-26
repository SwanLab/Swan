function CECMpostprocessCARTESIAN(xCECM,wCECM,zCECM,zDECM,wDECM,DATA,DATAapprox,xINI,A) 

if nargin == 0
    load('tmp.mat')
end

ndim = size(xCECM,2) ;

% ----------------------------------------------------------------------
%%% CHECKING THAT OUR RULE CAN INTEGRATE THE ORIGINAL SET OF FUNCTIONS 
% ----------------------------------------------------------------------
VSinv = [] ;
EVALUATE_GRADIENT = 0 ; 
DATA = DefaultField(DATA,'Evaluate_Fun_Gradients_via_FITTING',0) ; 
if DATA.Evaluate_Fun_Gradients_via_FITTING == 1
    DATA.PLOT_functions_SVD_decompositionFITTING = 0 ; 
    DATAFITTING = GetCoefficientesFitting(xINI,A',DATA) ; 
else
    DATAFITTING = [] ; 
end

[fCECM,~ ]= Evaluate_Basis_Grad_Analytic(xCECM,VSinv,DATA,EVALUATE_GRADIENT,DATAFITTING)  ;
EXACT_INTEGRAL = DATAapprox.ExactIntegral ; % Xf'*W ;
APPR_INTEGRAL  = fCECM*wCECM ;
errorMINE = norm(EXACT_INTEGRAL-APPR_INTEGRAL)./norm(EXACT_INTEGRAL)*100 ;
disp('-------------------------------------------------------------------------------------------')
disp(['Integration Error (original function) using OUR RULE WITH m = ',num2str(length(wCECM)),' points = ',num2str(errorMINE),' %']) ;
disp('-------------------------------------------------------------------------------------------')
%DATA =DefaultField(DATA,'NO_GAUSSIAN_RULE_COMPARISON',0) ; 
%if DATA.NO_GAUSSIAN_RULE_COMPARISON == 0
    % ----------------------------------------------------------
    % Comparison with  Standard GAUSS RULE (for polynomials)
    %----------------------------------------------------------
    [wGAUSS,xGAUSS,errorGAUSS] = StandardGaussRule(DATA,EXACT_INTEGRAL,size(xINI,2),DATAFITTING) ;
    DATAadd.wGAUSS = wGAUSS; DATAadd.xGAUSS = xGAUSS; DATAadd.errorGAUSS = errorGAUSS ;
   
%else
   % DATAadd.wGAUSS = [] ;  DATAadd.errorCECM = [] ; DATAadd.errorGAUSS = [] ;
%end
 DATAadd.errorCECM = errorMINE ; DATAadd.errorDECM = DATAapprox.ErrorApproxDECM ;
if  ndim ==3
    Plot3DpointsCECMgen(DATA,xINI(zDECM,:),xCECM,xINI,wCECM,wDECM,zCECM,DATAadd,DATAapprox,zDECM) ;
elseif ndim ==2
    Plot2DpointsCECMgen(DATA,xINI(zDECM,:),xCECM,xINI,wCECM,wDECM,zCECM,DATAadd,DATAapprox,zDECM) ;
elseif ndim == 1
    
    Plot1DpointsCECMgen(DATA,xINI(zDECM,:),xCECM,xINI,wCECM,wDECM,zCECM,DATAadd,DATAapprox,zDECM) ;
    
end




% % Plot Distorted domain
% if  DATAFUN.TYPE ==12
%     PlotDistortedDomainGauss3D(DATAOUT_irreg_gauss,DATAOUT_gauss,wGAUSS,xGAUSS,VSinv,DATALOC,xMAT,errorGAUSS,w_gauss_q)  ;
% end

