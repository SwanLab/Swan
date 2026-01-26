
%load('tmp1.mat')
DATAcoarsefine = [] ;
% INPUTS_LOC = DefaultField(INPUTS_LOC,'FACTOR_PERIODICITY',1) ;
% beta_p = INPUTS_LOC.FACTOR_PERIODICITY ;
a  = INPUTS_LOC.DISPLACEMENTS_RIGID_BODY ; % Rigid body amplitudes face 1


[Gb,dR,DOFr,DOFm] = FIXED_FACES_LINEAR_LATERAL_fun(a,DOMAINVAR,COOR,CONNECTb,TypeElementB,DATA) ;
