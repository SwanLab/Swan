function [setPoints,wRED,ERROR_GLO,DATAOUT] = DiscreteECMcable1D(BstRED_l,BasisTENSIONV,DATA,wSTs,DATAoffline)

 if nargin == 0
     load('tmp1.mat')
 end

% Basis matrix for internal forces 
% **********************************

 [Q,SNAPfint] = QbasisMatrixIntegrand_cable1D(BstRED_l,BasisTENSIONV,DATA,wSTs,DATAoffline) ; 

% Empirical cubature method
% -------------------------

%[setPoints,wRED,ERROR_GLO,DATAOUT]= EmpiricalCubatureMethod_orig(Q,[],wSTs,DATA_ECM)  ;

% CANDIDATE POINTS EXCLUDED FROM INFO IN FILE 
% -------------------------------------------

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



IntegralExact =SNAPfint'*wSTs  ; 
IntegrationError = IntegralExact - SNAPfint(setPoints,:)'*wRED ;
IntegrationError = norm(IntegrationError)/norm(IntegralExact); 
disp(['Actual integration error using DECM= ',num2str(IntegrationError*100),' % (prescribed tolerance fint =',num2str(DATAoffline.errorFINT*100), '%'])