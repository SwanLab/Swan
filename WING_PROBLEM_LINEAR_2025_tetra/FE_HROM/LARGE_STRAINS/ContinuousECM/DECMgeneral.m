function [ECMdata,HYPERREDUCED_VARIABLES,DATAOUTdecm] = DECMgeneral(A,U,wFE,xFE,DATA_ECM,MESH,DATA)


% List of elements to be excluded from the initial SET
DATA_ECM = DefaultField(DATA_ECM,'ListElementsToExclude',[]) ;
INDSEL = 1:length(wFE) ;
if ~isempty(DATA_ECM.ListElementsToExclude)
    ListGaussToExclude = small2large(DATA_ECM.ListElementsToExclude,MESH.ngausE) ;
    INDSEL  = setdiff(INDSEL,ListGaussToExclude) ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Discrete Empirical cubature method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DATA_ECM = DefaultField(DATA_ECM,'errorDECM',0) ;
DATA_DECM = [] ;
DATA_DECM.TOL = DATA_ECM.errorDECM ;
DATA_DECM.IND_POINTS_CANDIDATES = INDSEL ;
DATA_DECM.STORE_INFO_ITERATIONS = 1; 
[zDECM,wDECM,~,DATAOUTdecm]= DiscreteEmpiricalCubatureMethod(U',wFE,DATA_DECM)  ;
ECMdata.xDECM = xFE(zDECM,:);
ECMdata.wDECM = wDECM ;
% Determining the indices of the associated elements
setElements = large2smallREP(zDECM,MESH.ngausE) ;
disp('****************************+')
disp(['List of selected m = ',num2str(length(setElements)),' elements'])
disp(num2str(setElements'))
%clipboard('copy',num2str(setElements'));

HYPERREDUCED_VARIABLES.setPoints = zDECM ;  % SEt integration points
HYPERREDUCED_VARIABLES.setElements = setElements ;  % Set associated elements
HYPERREDUCED_VARIABLES.WdomRED = wDECM ;  % Set associated WEights
HYPERREDUCED_VARIABLES.PHI = bsxfun(@times,U,1./sqrt(wFE)) ;


% ------------------------------------------------
% Computing integration errors due to the DECM
% -----------------------------------------------

DATA = DefaultField(DATA,'CalculateErrorIntegralSnapshotMatrices',1) ; % = 0 ;

if DATA.CalculateErrorIntegralSnapshotMatrices == 1
    if isempty(DATA.ExactIntegral)
        error('Adapt the code so that it can use ErrorApproximateIntegral2 (more effective)')
        ErrorApproximateIntegral(A,U,DATA.ExactIntegral,wFE,DATA_ECM,zDECM,wDECM,DATA)  ;   % Before 25-Oct-2022
    else
        ErrorApproximateIntegral2(A,HYPERREDUCED_VARIABLES.PHI,wFE,DATA_ECM,zDECM,wDECM,DATA)  ;
    end
end
