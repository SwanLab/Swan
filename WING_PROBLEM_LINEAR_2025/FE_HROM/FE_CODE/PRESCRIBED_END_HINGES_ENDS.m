
a_A = INPUTS_LOC.DISPLACEMENTS_ENDS_RIGID_BODY{1} ; % Rigid body amplitudes face 1
a_B = INPUTS_LOC.DISPLACEMENTS_ENDS_RIGID_BODY{2} ; % Rigid body amplitudes face 2çç


[Gb,dR,DOFr,DOFm,AREA,R,DOFA,DOFB] = ...
    PRESCRIBED_END_HINGES_ENDSfun(a_A,a_B,DOMAINVAR,COOR,CONNECTb,TypeElementB,...
    DATA) ; 