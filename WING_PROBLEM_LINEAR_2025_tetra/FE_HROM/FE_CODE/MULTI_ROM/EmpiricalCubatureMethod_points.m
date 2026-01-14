function [HYPERREDUCED_VARIABLES,MSG,Wdom] = EmpiricalCubatureMethod_points(DATAIN,DATAROM,BASES,DATA_REFMESH)

% Compute basis matrix for internal forces. We need the basis matrix for the stresses, the associated singular values
% as well as the reduced Bmatrix
% -----------------------------------------------------------------------
% Bases for displacements
% BasisU = DATAROM.BasisUdef ;
% Reduced B matrix
if nargin == 0
    load('tmp2.mat')
    %  DATAIN.CUBATURE.AlignStressesWithStrains = 1;
end
disp('****************************+')
disp('Empirical Cubature Method ')
disp('*************************************')
load(DATAIN.NAME_WS_MODES,'CgloDOM','Wdom','Bdom')



DATAIN = DefaultField(DATAIN,'CUBATURE',[]) ;
DATAIN.CUBATURE = DefaultField(DATAIN.CUBATURE,'IncludeSingularValues_displacements',0) ; % Include Singular Value displacement
MSG = {} ;

BdomRED = Bdom*DATAROM.BasisUdef ;

% Methods for computing internal forces snapshots
DATAIN.CUBATURE = DefaultField(DATAIN.CUBATURE,'INTERNAL_FORCES_WITHOUT_STRESS_MODES',0);
if DATAIN.CUBATURE.INTERNAL_FORCES_WITHOUT_STRESS_MODES  == 1
    % Method implemented in 26-Jan-2020. Internal forces computed directly,
    % without the need for stress modes
    [BasisF, SingVal_F,MSG] = InternalFocesECM_WithoutStressModes(BASES,DATAIN,MSG,BdomRED,Wdom) ;
    
    
else
    % Classic method (from stress modes)
    [BasisF, SingVal_F,MSG] = InternalFocesECM_fromStressModes(BASES,DATAIN,MSG,BdomRED,Wdom,DATAROM,Bdom) ;
end


%     % Empirical Cubature Method
DATA_CUBATURE = []  ;
DATAIN.CUBATURE = DefaultField(DATAIN.CUBATURE,'IncludeSingularValuesF',0) ;
DATA_CUBATURE.IncludeSingularValuesF = DATAIN.CUBATURE.IncludeSingularValuesF ;

DATAIN.CUBATURE = DefaultField(DATAIN.CUBATURE,'TOLERANCE_EmpiricalCubature',0) ;

DATA_CUBATURE.TOL = DATAIN.CUBATURE.TOLERANCE_EmpiricalCubature ;

% The  number of integrations points cannot be lower than the numbe

[setPoints,WdomRED]= EmpiricalCubatureMethod(BasisF,SingVal_F,Wdom,DATA_CUBATURE) ;
[setPoints ixxx]= sort(setPoints) ;
WdomRED = WdomRED(ixxx) ;
% Determining the indices of the associated elements
ngausTOTAL = length(Wdom);
nstrain = size(BdomRED,1)/length(Wdom);
ngaus = ngausTOTAL/size(DATA_REFMESH.CN,1) ;
setElements = large2smallREP(setPoints,ngaus) ;
setElements_WITH_REPET = large2smallINCLUDEREP(setPoints,ngaus) ;

disp('****************************+')
disp(['List of selected m = ',num2str(length(setElements)),' elements'])
disp(num2str(setElements'))
clipboard('copy',num2str(setElements'));

MSG{end+1} = ['List of selected m = ',num2str(length(setElements)),' elements'] ;
MSG{end+1} = num2str(setElements') ;
% Computation reduced-order stiffness matrix , as well as r
DATAloc = [] ; 
[Kstiff_hyper,Celas_Bdom,CelasBdomALL] = StiffReduced(setPoints, nstrain,CgloDOM,DATAIN,BdomRED,WdomRED,Wdom,DATAloc) ;




% Comparison with the one computed used the whole set of Gauss points
Kstiff_all = DATAROM.KdomRED ;
DIFFERENCE = norm((Kstiff_hyper-Kstiff_all),'fro')/norm(Kstiff_all,'fro')*100 ;
MSGloc1 = ['Difference between HROM and FE Stiff Matrix = ',num2str(DIFFERENCE),' %'] ;
disp(MSGloc1)
MSG{end+1} = [MSGloc1,'(per cent)'] ;
disp('*************************************************************************************+')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HYPERREDUCED_VARIABLES.setPoints = setPoints ;  % SEt integration points
HYPERREDUCED_VARIABLES.setElements = setElements ;  % Set associated elements
HYPERREDUCED_VARIABLES.setElements_WITH_REPET = setElements_WITH_REPET ;  % Set associated elements

HYPERREDUCED_VARIABLES.WdomRED = WdomRED ;  % Set associated WEights
HYPERREDUCED_VARIABLES.Celas_Bdom = Celas_Bdom ;  % Product Elastic Matrix times  B matrix
HYPERREDUCED_VARIABLES.nstrain = nstrain ;
HYPERREDUCED_VARIABLES.CelasBdomALL = CelasBdomALL ;


end
