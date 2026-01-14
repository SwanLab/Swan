
%load('tmp1.mat')
DATAcoarsefine = [] ;
% INPUTS_LOC = DefaultField(INPUTS_LOC,'FACTOR_PERIODICITY',1) ;
% beta_p = INPUTS_LOC.FACTOR_PERIODICITY ;
a  = INPUTS_LOC.DISPLACEMENTS_RIGID_BODY ; % Rigid body amplitudes face 1

%INPUTS_LOC = DefaultField(INPUTS_LOC,'NameFileMeshLOC_coarse',INPUTS_LOC.NameFileMesh) ;
%

%OKK = strcmp(INPUTS_LOC.NameFileMeshLOC_coarse,INPUTS_LOC.NameFileMesh) ;
%if OKK ==1
    [Gb,dR,DOFr,DOFm,NODES_ENTITIES,ANG_ROTATION_TOTAL] = FIXED_FACES_RVES_fun(a,DOMAINVAR,COOR,CONNECTb,TypeElementB,DATA) ;
     
%else
    
%    error('Option not implemented')
    
%     [Gb,dR,DOFr,DOFm,AREA,R,DOFA,DOFB,DATAcoarsefine] = ...
%         PERIODIC_BEAMS_COARSE_FINE(a_A,a_B,DOMAINVAR,COOR,CONNECTb,...
%         TypeElementB,INPUTS_LOC.NameFileMeshLOC_coarse,DATA) ;
    
%end


%DOMAINVAR.RigidBodyModes =R ;


