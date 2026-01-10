
%load('tmp1.mat')

% INPUTS_LOC = DefaultField(INPUTS_LOC,'FACTOR_PERIODICITY',1) ; 
% beta_p = INPUTS_LOC.FACTOR_PERIODICITY ; 
a_A = INPUTS_LOC.DISPLACEMENTS_ENDS_RIGID_BODY{1}(:) ; % Rigid body amplitudes face 1
a_B = INPUTS_LOC.DISPLACEMENTS_ENDS_RIGID_BODY{2}(:) ; % Rigid body amplitudes face 2

b_A = INPUTS_LOC.GENERALIZED_FORCES_ENDS_BEAM{1}(:) ; % Forces amplitudes face 1
b_B = INPUTS_LOC.GENERALIZED_FORCES_ENDS_BEAM{2}(:) ; % Forces amplitudes face 2
if iscell(a_A)
    a_A = cell2mat(a_A) ; 
end
if iscell(a_B)
    a_B = cell2mat(a_B) ; 
end
INPUTS_LOC = DefaultField(INPUTS_LOC,'NameFileMeshLOC_coarse',INPUTS_LOC.NameFileMesh) ; 

%OKK = strcmp(INPUTS_LOC.NameFileMeshLOC_coarse,INPUTS_LOC.NameFileMesh) ; 
%if OKK ==1
[Gb,dR,DOFr,DOFm,AREA,R,Fpnt_Dirichlet] = ...
    BOUNDARY_COND_SIMPLE_BENDING_fun(a_A,a_B,DOMAINVAR,COOR,CONNECTb,TypeElementB,DATA,b_A,b_B) ; 





% else
%  [Gb,dR,DOFr,DOFm,AREA,R] = ...
%      PERIODIC_BEAMS_COARSE_FINE(a_A,a_B,DOMAINVAR,COOR,CONNECTb,...
%      TypeElementB,INPUTS_LOC.NameFileMeshLOC_coarse,DATA) ;  ; 
%    
% end

DOMAINVAR.RigidBodyModes =R ; 

    
    