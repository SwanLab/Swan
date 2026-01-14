function  ErrorCalcLocal2023(DATA,CECMoutput,DATA_AUX,A,DATAloc)

if nargin == 0
    load('tmp1.mat')
end


% Approximated integral
DATALOC=[] ;

DATA_AUX.VAR_SMOOTH_FE.BasisIntegrand = A ;
DATA_AUX.VAR_SMOOTH_FE.ExactIntegral = DATAloc.ExactIntegral ;

[Acecm,~,POLYINFO]=     EVALBASIS(CECMoutput.xCECM,DATALOC,DATA_AUX.VAR_SMOOTH_FE,DATA_AUX.POLYINFO)  ;
IntApproxCECM = Acecm'*CECMoutput.wCECM ;



ErrorApprox = norm(DATAloc.ExactIntegral-IntApproxCECM)/norm(DATAloc.ExactIntegral)*100 ;
disp(['Approximation ERROR CECM original functions  (m= ',num2str(length(CECMoutput.wCECM)),') = ',num2str(ErrorApprox),' %']) ;


DATA = DefaultField(DATA,'xGAUSSIAN',[]) ;

if  ~isempty(DATA.xGAUSSIAN)
    DATA_AUX.POLYINFO.setElements = [] ; 
    [Agauss,~,POLYINFO]=     EVALBASIS(DATA.xGAUSSIAN,DATALOC,DATA_AUX.VAR_SMOOTH_FE,DATA_AUX.POLYINFO)  ;
    
    if size(Agauss,1) == length(DATA.weigthsGAUSSIAN)
    IntApproxGAUSS= Agauss'*DATA.weigthsGAUSSIAN ;
   
    
    %         vGAUSS = sum(DATA.xGAUSSIAN) ;
    %         vCECM = sum(CECMoutput.wCECM) ;
    % %         if abs(vGAUSS-vCECM)/abs(vCECM) > 1e-2
    %             warning('Gauss integration might be ill-calculated')
    %         end
    
    ErrorApprox = norm(DATAloc.ExactIntegral-IntApproxGAUSS)/norm(DATAloc.ExactIntegral)*100 ;
    disp(['Approximation ERROR GAUSS original functions  (m= ',num2str(length(DATA.weigthsGAUSSIAN)),') =',num2str(ErrorApprox),' %'])  ;
    
    
    disp(['COORDINATE and WEIGHTS GAUSSIAN POINTS'])
    disp('----------------------------------------------')
    x = DATA.xGAUSSIAN(:,1) ;
    y = DATA.xGAUSSIAN(:,2) ;
    w = DATA.weigthsGAUSSIAN  ;
    table(x,y,w)
    
    else
       disp('Gaussian rule cannot be computed (empty domain)') 
        
    end
    
    
end
disp('----------------------------------------------')


disp(['COORDINATE and WEIGHTS CECM POINTS'])
disp('----------------------------------------------')
x = CECMoutput.xCECM(:,1) ;
y = CECMoutput.xCECM(:,2) ;
w = CECMoutput.wCECM ;
table(x,y,w)
