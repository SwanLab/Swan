function [HYPERREDUCED_VARIABLES,MSG ] =...
    ReducedSetIntegrationPoints(DATAIN,DATAROM,BASES,DATA_REFMESH,MSG)

if nargin == 0
    load('tmp2.mat')
    %  DATAIN.CUBATURE.USE_INTERPOLATION_POINTS_STRESSES = 2;
end


DATAIN = DefaultField(DATAIN,'ISNONLINEAR',0) ;
DATAIN = DefaultField(DATAIN,'CUBATURE',[]) ;

% if DATAIN.ISNONLINEAR == 0
DATAIN.CUBATURE = DefaultField(DATAIN.CUBATURE,'USE_INTERPOLATION_POINTS_STRESSES',0) ;
% else
%     DATAIN.CUBATURE = DefaultField(DATAIN.CUBATURE,'USE_INTERPOLATION_POINTS_STRESSES',0) ;
%
% end


MSGloc = [] ;
if  DATAIN.CUBATURE.USE_INTERPOLATION_POINTS_STRESSES == 0
    [HYPERREDUCED_VARIABLES,MSGloc,Wdom] = EmpiricalCubatureMethod_points(DATAIN,DATAROM,BASES,DATA_REFMESH) ;
elseif  DATAIN.CUBATURE.USE_INTERPOLATION_POINTS_STRESSES == 1
    HYPERREDUCED_VARIABLES = StressPoints(DATAIN,DATAROM,BASES,DATA_REFMESH) ;
elseif DATAIN.CUBATURE.USE_INTERPOLATION_POINTS_STRESSES == 2
    HYPERREDUCED_VARIABLES = StressPointsCubature(DATAIN,DATAROM,BASES,DATA_REFMESH) ;
    
else
    error('Option not implemented')
end

for iii = 1:length(MSGloc)
    MSG{end+1} = MSGloc{iii} ;
end

% Reconstruction matrix
% ---------------------
DATAIN = DefaultField(DATAIN,'COMPUTE_RECONSTRUCTION_MATRIX',0) ;
if  DATAIN.COMPUTE_RECONSTRUCTION_MATRIX  == 1 || DATAIN.ISNONLINEAR==1
    BasisS = BASES.STRESSES.U  ;
    setGauss = small2large(HYPERREDUCED_VARIABLES.setPoints,HYPERREDUCED_VARIABLES.nstrain) ;
    coeff = (BasisS(setGauss,:)'*BasisS(setGauss,:))\BasisS(setGauss,:)' ;
    HYPERREDUCED_VARIABLES.ReconsStresses = BasisS*coeff ;
    
    % HOMOGENIZATION OPERATOR  (19-JAN-2020, for the paper)
    % --------------------------------------------------------------------------------------------
    load(DATAIN.NAME_WS_MODES,'Wdom')
    nstrain = size(BasisS,1)/length(Wdom); 
    HOMOGENIZATION_OPERATOR = zeros(nstrain,size(HYPERREDUCED_VARIABLES.ReconsStresses,2)) ; 

    for istrain = 1:nstrain 
        COMP = istrain:nstrain:size(HYPERREDUCED_VARIABLES.ReconsStresses,1) ; 
        HOMOGENIZATION_OPERATOR(istrain,:) = Wdom'*HYPERREDUCED_VARIABLES.ReconsStresses(COMP,:) ; 
    end
   HYPERREDUCED_VARIABLES.HOMOGENIZATION_OPERATOR = HOMOGENIZATION_OPERATOR ; 
   % ------------------------------------------------------------------------------------------- 
    
end

% Printing reduced set of elements 
%--------------------------------------
MSG{end+1} = '-----------------------------------------' ; 
MSG{end+1} = 'Position ECM points (ans weights )' ; 
MSG{end+1} = '-----------------------------------------' ; 

DATAIN = PrintingGid_ECMpoints(DATAIN,DATA_REFMESH,HYPERREDUCED_VARIABLES) ;

for iprint =  1:length(DATAIN.MSGPRINT)
    MSG{end+1} = DATAIN.MSGPRINT{iprint} ; 
end

% -------------------------------------
MSG{end+1} = '-----------------------------------------' ; 



end



function HYPERREDUCED_VARIABLES = StressPoints(DATAIN,DATAROM,BASES,DATA_REFMESH)

% Compute basis matrix for internal forces. We need the basis matrix for the stresses, the associated singular values
% as well as the reduced Bmatrix
% -----------------------------------------------------------------------
% Bases for displacements
% BasisU = DATAROM.BasisUdef ;
% Reduced B matrix
disp('Computing stress points')
disp('************************************++')
load(DATAIN.NAME_WS_MODES,'CgloDOM','Wdom','Bdom')
BdomRED = Bdom*DATAROM.BasisUdef ;
BasisS = BASES.STRESSES.U ;
ndim = size(DATA_REFMESH.COOR,2) ;
nstrain = 6 ;
if ndim == 2
    nstrain = 4;
end
%
npoints = size(BasisS,1)/nstrain ;
BasisScomp = zeros(npoints,size(BasisS,2)*nstrain) ;

% Split into components
for istrain = 1:nstrain
    
    IND = istrain:nstrain:size(BasisS,1) ;
    COLS = istrain:nstrain:size(BasisScomp,2) ;
    BasisScomp(:,COLS) = BasisS(IND,:) ;
    
end

%% Singular value decomposition
[U,S] = RSVDT(BasisScomp) ;
% DEIM
setPoints = DEIMfun(U) ;
[setPoints ixxx]= sort(setPoints) ;
% Determining the indices of the associated elements
ngaus = size(BasisScomp,1)/size(DATA_REFMESH.CN,1) ;
setElements = large2small(setPoints,ngaus) ;
disp('****************************+')
disp(['List of selected m = ',num2str(length(setElements)),' elements'])
disp(num2str(setElements'))
clipboard('copy',num2str(setElements'));
% Computation  of the reduced stiffness using this set of points
% ---------------------------------------------------------------
%  Set of entries associated with these set of points
setGauss = small2large(setPoints,nstrain) ;
% HyperReduced B matrix
BdomRED_hyper = BdomRED(setGauss,:) ;
% HyperReduced elasticity matrix  ---including full-order weights
wCelas_hyper = CgloDOM(setGauss,setGauss) ;
% Hyperreduced matrix multiplied by weights
% PRoduct C*B
wCelas_Bdom = wCelas_hyper*BdomRED_hyper ;
CelasBdom = zeros(size(wCelas_Bdom)) ;
for istrain = 1:nstrain
    INDS = istrain:nstrain:length(setGauss) ;
    Celas_Bdom(INDS,:) = bsxfun(@times,wCelas_Bdom(INDS,:),1./Wdom(setPoints)) ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HYPERREDUCED_VARIABLES.setPoints = setPoints ;  % SEt integration points
HYPERREDUCED_VARIABLES.setElements = setElements ;  % Set associated elements
HYPERREDUCED_VARIABLES.WdomRED = [] ;  % Set associated WEights
HYPERREDUCED_VARIABLES.Celas_Bdom = Celas_Bdom ;  % Product Elastic Matrix times  B matrix
HYPERREDUCED_VARIABLES.nstrain = nstrain ;

end



function HYPERREDUCED_VARIABLES = StressPointsCubature(DATAIN,DATAROM,BASES,DATA_REFMESH)

% Compute basis matrix for internal forces. We need the basis matrix for the stresses, the associated singular values
% as well as the reduced Bmatrix
% -----------------------------------------------------------------------
% Bases for displacements
% BasisU = DATAROM.BasisUdef ;
% Reduced B matrix

% Cubature points
% ----------------
HYPERREDUCED_VARIABLES = EmpiricalCubatureMethod_points(DATAIN,DATAROM,BASES,DATA_REFMESH) ;
setPointsCub = HYPERREDUCED_VARIABLES.setPoints ;
setElementsCub = HYPERREDUCED_VARIABLES.setElements ;
% Stress points
% ..............
disp('Computing stress points')
disp('************************************++')
load(DATAIN.NAME_WS_MODES,'CgloDOM','Wdom','Bdom')
BdomRED = Bdom*DATAROM.BasisUdef ;
BasisS = BASES.STRESSES.U ;
ndim  = size(DATA_REFMESH.COOR,2) ;
if ndim == 3
    nstrain = 6 ;
else
    nstrain = 4;
end
%
npoints = size(BasisS,1)/nstrain ;
BasisScomp = zeros(npoints,size(BasisS,2)*nstrain) ;

% Split into components
for istrain = 1:nstrain
    
    IND = istrain:nstrain:size(BasisS,1) ;
    COLS = istrain:nstrain:size(BasisScomp,2) ;
    BasisScomp(:,COLS) = BasisS(IND,:) ;
    
end

%% Singular value decomposition
[U,S] = RSVDT(BasisScomp) ;
% DEIM
setPoints = DEIMfun(U) ;

setPoints = unique([setPoints;setPointsCub]) ;

[setPoints ixxx]= sort(setPoints) ;
% Determining the indices of the associated elements
ngaus = size(BasisScomp,1)/size(DATA_REFMESH.CN,1) ;
setElements = large2small(setPoints,ngaus) ;
disp('****************************+')
disp(['List of selected m = ',num2str(length(setElements)),' elements'])
disp(num2str(setElements'))
clipboard('copy',num2str(setElements'));
% Computation  of the reduced stiffness using this set of points
% ---------------------------------------------------------------
%  Set of entries associated with these set of points
setGauss = small2large(setPoints,nstrain) ;
% HyperReduced B matrix
BdomRED_hyper = BdomRED(setGauss,:) ;

% HyperReduced elasticity matrix  ---including full-order weights
wCelas_hyper = CgloDOM(setGauss,setGauss) ;
% Hyperreduced matrix multiplied by weights
% PRoduct C*B
wCelas_Bdom = wCelas_hyper*BdomRED_hyper ;
CelasBdom = zeros(size(wCelas_Bdom)) ;
for istrain = 1:nstrain
    INDS = istrain:nstrain:length(setGauss) ;
    Celas_Bdom(INDS,:) = bsxfun(@times,wCelas_Bdom(INDS,:),1./Wdom(setPoints)) ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HYPERREDUCED_VARIABLES.setPoints = setPoints ;  % SEt integration points
HYPERREDUCED_VARIABLES.setElements = setElements ;  % Set associated elements
HYPERREDUCED_VARIABLES.WdomRED = [] ;  % Set associated WEights
HYPERREDUCED_VARIABLES.Celas_Bdom = Celas_Bdom ;  % Product Elastic Matrix times  B matrix
HYPERREDUCED_VARIABLES.setPointsCub = setPointsCub ;  % SEt integration points
HYPERREDUCED_VARIABLES.setElementsCub = setElementsCub ;  % Set associated elements
HYPERREDUCED_VARIABLES.nstrain = nstrain ;

end




