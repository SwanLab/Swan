function [setPoints,wRED,ERROR_GLO,DATAOUT] = DiscreteECM(BstRED_l,BasisPone,DATA,wSTs,DATAoffline)

 

% Basis matrix for internal forces 
% **********************************
TIMEq = tic ; 
 [Q,SNAPfint] = QbasisMatrixIntegrand(BstRED_l,BasisPone,DATA,wSTs,DATAoffline) ; 
DATAOUT.TIMEq = toc(TIMEq) ; 

    
    disp(['Number of internal force modes = ',num2str(size(Q,2))])


disp(['Time to assembly orthogonal matrix of the integrand =',num2str(DATAOUT.TIMEq)])



% Empirical cubature method
% -------------------------

%[setPoints,wRED,ERROR_GLO,DATAOUT]= EmpiricalCubatureMethod_orig(Q,[],wSTs,DATA_ECM)  ;

% CANDIDATE POINTS EXCLUDED FROM INFO IN FILE 
% -------------------------------------------
TIMEdecm = tic ; 

DATAoffline = DefaultField(DATAoffline,'ListElementsExclude_fromGID',[]) ; 
if ~isempty(DATAoffline.ListElementsExclude_fromGID)
     ListElementsToExclude = load(DATAoffline.ListElementsExclude_fromGID) ; 
    ListElementsToExclude = ListElementsToExclude(:,1) ; 
    ngausELEM = DATA.MESH.ngaus_STRESS; 
    ListGaussToExclude = small2large(ListElementsToExclude,ngausELEM) ; 
    INDSEL  = setdiff(1:size(Q,1),ListGaussToExclude) ; 
else
    INDSEL  = 1:size(Q,1) ; 
end

DATA_ECM.IND_POINTS_CANDIDATES = INDSEL ;
DATA_ECM.TOL = DATAoffline.errorECM ;
DATAoffline = DefaultField(DATAoffline,'USE_SELECTIVE_DECM',1) ; % = 1;
if DATAoffline.USE_SELECTIVE_DECM == 0
    % Version before 3-DEc-2021
    [setPoints,wRED,ERROR_GLO,DATAOUT]= EmpiricalCubatureMethod_orig(Q,[],wSTs,DATA_ECM)  ;
else
    % New version (after 3-DEc-2021)
    [setPoints,wRED,ERROR_GLO,DATAOUT]= DiscreteEmpiricalCubatureMethod(Q',wSTs,DATA_ECM)  ;
    
end

DATAOUT.TIMEdecm = toc(TIMEdecm) ; 

disp(['Time to select the integration points =',num2str(DATAOUT.TIMEdecm)])




IntegralExact =SNAPfint'*wSTs  ; 
IntegrationError = IntegralExact - SNAPfint(setPoints,:)'*wRED ;
IntegrationError = norm(IntegrationError)/norm(IntegralExact); 
disp(['Actual integration error using DECM= ',num2str(IntegrationError*100),' % (prescribed tolerance fint =',num2str(DATAoffline.errorFINT*100), '%'])