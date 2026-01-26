function  [xCECM,wCECM,zCECM,zDECM,wDECM,DATAapprox] =  ContinuousEmpiricalCubatureMethod_givenDECM(PHI,W,xINI,zDECM,wDECM,DATA)
% Adaptation ContinuousEmpiricalCubatureMethod.m
% JAHO, 16-Dec-2021 
% --------------------------------------------------------
if nargin == 0
    load('tmp.mat')
end

% STEP 1) 
% -------
% SVD of A*sqrt(W)
% -----------------------------
%DATA = DefaultField(DATA,'TOLsvdXf',0); % Tolerance SVD 
%DATA = DefaultField(DATA,'RELOAD_SNAPSHOT_MATRIX',-1); % Load previous results 
%if DATA.RELOAD_SNAPSHOT_MATRIX == 0 || DATA.RELOAD_SNAPSHOT_MATRIX == -1
%     [GbarT,S_A,V_A,INTexac,DATA] = ObtainBasisMatrix(A',W,DATA) ; 
%    if DATA.RELOAD_SNAPSHOT_MATRIX == 0
%        save(DATA.NAMEWS,'GbarT','S_A','V_A','INTexac','-append')
%    end
%else
%    load(DATA.NAMEWS,'GbarT','S_A','V_A','INTexac')
%end

GbarT = bsxfun(@times,PHI,sqrt(W)) ; 
S_A = [] ; 
V_A = [] ; 
INTexac = PHI'*W ; 

% STEP 2) Discrete Empirical Cubature Method (initial set of integration points )
% ------------------------------------------
% Points included in the initial search set   
% DATA = DefaultField(DATA,'ECM_POINTS_INCLUDE',[]) ; 
% DATA = DefaultField(DATA,'TOL_DECM',0) ; 
% DATAecm.IND_POINTS_CANDIDATES = DATA.ECM_POINTS_INCLUDE ; 
% DATAecm.TOL = DATA.TOL_DECM ; 
 %DATAecm.SingularValuesSnapshotMatrix =S_A ; 


%[zDECM,wDECM,~,~]= DiscreteEmpiricalCubatureMethod(GbarT',W,DATAecm)  ;
% Approximation Error original matrix using the DECM points  
ErrorApproxDECM = ShowApproximationERROR_DECM(PHI,zDECM,wDECM,INTexac)  ; 
  
% STEP 3) Removal algorithm  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[xCECM,wCECM,zCECM,DATAapprox] = RemoveRecursivePointsCART(GbarT',W,xINI,zDECM,wDECM,DATA,S_A,V_A) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 

DATAapprox.ExactIntegral = INTexac ;
DATAapprox.ErrorApproxDECM = ErrorApproxDECM ;

