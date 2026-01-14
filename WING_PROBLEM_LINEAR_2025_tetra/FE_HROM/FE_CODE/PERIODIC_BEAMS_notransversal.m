%load('tmp1.mat')
a_A = INPUTS_LOC.DISPLACEMENTS_ENDS_RIGID_BODY{1}(:) ; % Rigid body amplitudes face 1
a_B = INPUTS_LOC.DISPLACEMENTS_ENDS_RIGID_BODY{2}(:) ; % Rigid body amplitudes face 2
if iscell(a_A)
    a_A = cell2mat(a_A) ; 
end
if iscell(a_B)
    a_B = cell2mat(a_B) ; 
end
%da = a_A-a_B ;
[Gb,dR,DOFr,DOFm]  = PeriodicBCs_unrestrictedTRANSVERSAL(DOMAINVAR,COOR,CONNECTb,TypeElementB,a_A,a_B) ; 