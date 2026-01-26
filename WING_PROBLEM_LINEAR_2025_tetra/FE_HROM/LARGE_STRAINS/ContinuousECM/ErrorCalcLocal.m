function  ErrorCalcLocal(DATA,xLIM,CECMoutput) 


% Approximated integral
    Acecm = feval(DATA.Integrand.NameFunctionGenerate,xLIM,CECMoutput.xCECM,DATA.Integrand) ;
    IntApproxCECM = Acecm'*CECMoutput.wCECM ;
    
    ErrorApprox = norm(CECMoutput.INTexac-IntApproxCECM)/norm(CECMoutput.INTexac)*100 ;
    disp(['Approximation ERROR CECM original functions  (m= ',num2str(length(CECMoutput.wCECM)),') = ',num2str(ErrorApprox),' %']) ;
    
    
    DATA = DefaultField(DATA,'xGAUSS',[]) ;
    
    if  ~isempty(DATA.xGAUSS)
        Agauss = feval(DATA.Integrand.NameFunctionGenerate,xLIM,DATA.xGAUSS,DATA.Integrand) ;
        IntApproxGAUSS= Agauss'*DATA.wGAUSS ;
        
        vGAUSS = sum(DATA.wGAUSS) ;
        vCECM = sum(CECMoutput.wCECM) ;
        if abs(vGAUSS-vCECM)/abs(vCECM) > 1e-2
            warning('Gauss integration might be ill-calculated')
        end
        
        ErrorApprox = norm(CECMoutput.INTexac-IntApproxGAUSS)/norm(CECMoutput.INTexac)*100 ;
        disp(['Approximation ERROR GAUSS original functions  (m= ',num2str(length(DATA.wGAUSS)),') =',num2str(ErrorApprox),' %'])  ;
    end
    