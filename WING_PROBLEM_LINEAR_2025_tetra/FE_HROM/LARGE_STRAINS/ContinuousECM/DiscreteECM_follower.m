function [setPoints,wRED,ERROR_GLO,DATAOUT] = DiscreteECM_follower(NstRED_l,Basis_tPRESS,DATA,wSTs,DATAoffline)

if nargin == 0
    load('tmp.mat')
end
 

% Basis matrix for internal forces 
% **********************************

 Q = QbasisMatrixIntegrand_follower(NstRED_l,Basis_tPRESS,DATA,wSTs,DATAoffline) ; 

% Empirical cubature method
% -------------------------
DATA_ECM = [] ;
DATA_ECM.TOL = DATAoffline.errorECM ;
DATAoffline = DefaultField(DATAoffline,'TOLFilterCandidatePoints_ECM',1e-14) ;
DATA_ECM.TOLFilterCandidatePoints = DATAoffline.TOLFilterCandidatePoints_ECM ; 

[setPoints,wRED,ERROR_GLO,DATAOUT]= EmpiricalCubatureMethod_orig(Q,[],wSTs,DATA_ECM)  ;