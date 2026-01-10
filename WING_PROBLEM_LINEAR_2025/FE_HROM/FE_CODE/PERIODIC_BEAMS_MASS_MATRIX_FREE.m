% 
% %load('tmp1.mat')
% DATAcoarsefine = [] ; 
% % INPUTS_LOC = DefaultField(INPUTS_LOC,'FACTOR_PERIODICITY',1) ; 
% % beta_p = INPUTS_LOC.FACTOR_PERIODICITY ; 
% a_A = INPUTS_LOC.DISPLACEMENTS_ENDS_RIGID_BODY{1}(:) ; % Rigid body amplitudes face 1
% a_B = INPUTS_LOC.DISPLACEMENTS_ENDS_RIGID_BODY{2}(:) ; % Rigid body amplitudes face 2
% if iscell(a_A)
%     a_A = cell2mat(a_A) ; 
% end
% if iscell(a_B)
%     a_B = cell2mat(a_B) ; 
% end
% INPUTS_LOC = DefaultField(INPUTS_LOC,'NameFileMeshLOC_coarse',INPUTS_LOC.NameFileMesh) ; 
% 
% 
%  
%      [Gb,dR,DOFr,DOFm,AREA,R,DOFA,DOFB] = ...
%             PERIODIC_BEAMS_FREE(a_A,a_B,DOMAINVAR,COOR,CONNECTb,...
%             TypeElementB,INPUTS_LOC.NameFileMeshLOC_coarse,DATA) ;  
%  
% 
% DOMAINVAR.RigidBodyModes =R ; 
% 
%     
%     