
%load('tmp1.mat')

% INPUTS_LOC = DefaultField(INPUTS_LOC,'FACTOR_PERIODICITY',1) ; 
% beta_p = INPUTS_LOC.FACTOR_PERIODICITY ; 
a_A = INPUTS_LOC.DISPLACEMENTS_ENDS_RIGID_BODY{1}(:) ; % Rigid body amplitudes face 1
a_B = INPUTS_LOC.DISPLACEMENTS_ENDS_RIGID_BODY{2}(:) ; % Rigid body amplitudes face 2
if iscell(a_A)
    a_A = cell2mat(a_A) ; 
end
if iscell(a_B)
    a_B = cell2mat(a_B) ; 
end
INPUTS_LOC = DefaultField(INPUTS_LOC,'NameFileMeshLOC_coarse',INPUTS_LOC.NameFileMesh) ; 

%OKK = strcmp(INPUTS_LOC.NameFileMeshLOC_coarse,INPUTS_LOC.NameFileMesh) ; 
%if OKK ==1
[Gb,dR,DOFr,DOFm,AREA,R] = PRESCRIBED_FLUCTUATIONS_bending_fun(a_A,a_B,DOMAINVAR,COOR,CONNECTb,TypeElementB,DATA) ; 
% else
%  [Gb,dR,DOFr,DOFm,AREA,R] = ...
%      PERIODIC_BEAMS_COARSE_FINE(a_A,a_B,DOMAINVAR,COOR,CONNECTb,...
%      TypeElementB,INPUTS_LOC.NameFileMeshLOC_coarse,DATA) ;  ; 
%    
% end

DOMAINVAR.RigidBodyModes =R ; 

    
    