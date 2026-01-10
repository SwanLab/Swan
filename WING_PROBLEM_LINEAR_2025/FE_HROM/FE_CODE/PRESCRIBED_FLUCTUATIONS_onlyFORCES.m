% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/105_BackPeriodic/05_ONLY_FORCES.mlx
%load('tmp1.mat')
DATAcoarsefine = [] ;
% INPUTS_LOC = DefaultField(INPUTS_LOC,'FACTOR_PERIODICITY',1) ;
% beta_p = INPUTS_LOC.FACTOR_PERIODICITY ;
%a_A = INPUTS_LOC.DISPLACEMENTS_ENDS_RIGID_BODY{1}(:) ; % Rigid body
%amplitudes face 1  % NOT USED
% a_B = INPUTS_LOC.DISPLACEMENTS_ENDS_RIGID_BODY{2}(:) ; % Rigid body amplitudes face 2



U_b =  INPUTS_LOC.Ufluc_b ; 
%U_n =  INPUTS_LOC.Ufluc_n ; 

gammaB = INPUTS_LOC.gammaB ; 



% 
% D_amplitude_fluctuations =INPUTS_LOC.D_amplitude_fluctuations  ; 
% if iscell(a_A)
%     a_A = cell2mat(a_A) ;
% end
% if iscell(a_B)
%     a_B = cell2mat(a_B) ;
% end
INPUTS_LOC = DefaultField(INPUTS_LOC,'NameFileMeshLOC_coarse',INPUTS_LOC.NameFileMesh) ;


%%% CURVED ELEMENTS (rotation matrix FACE2 of trailing domain)
% -----------------
DATA = DefaultField(DATA,'angROTATION_FACE',[]) ;
DOMS = [1,size(CONNECTb,1)] ;
iface =2 ;
ROTMATRIX = LocalRotMatrix(iface,DATA,DOMS,ndim) ;


%if    ~isstruct(INPUTS_LOC.NameFileMeshLOC_coarse)
OKK = strcmp(INPUTS_LOC.NameFileMeshLOC_coarse,INPUTS_LOC.NameFileMesh) ;
if OKK ==1
    % No coarse mesh
    % ---------------
    if isempty(ROTMATRIX)
        % Straight domains
        [Gb,dR,DOFr,DOFm,AREA,R,DOFA,DOFB] = PRESCRIBED_FLUCTUATIONS_withFORCESonly(U_b,DOMAINVAR,COOR,CONNECTb,TypeElementB,...
            DATA,gammaB) ;
    else
        % Curved domains
        error('Option not implemented yet')
   %     [Gb,dR,DOFr,DOFm,AREA,R,DOFA,DOFB] = PERIODIC_BEAMS_MASS_MATRIXcurved(a_A,a_B,DOMAINVAR,COOR,CONNECTb,TypeElementB,DATA,ROTMATRIX) ;
    end
else
    % Fine and coarse mesh approach
    if isempty(ROTMATRIX)
        error('Option not implemented yet')
%         [Gb,dR,DOFr,DOFm,AREA,R,DOFA,DOFB,DATAcoarsefine] = ...
%             PERIODIC_BEAMS_COARSE_FINE(a_A,a_B,DOMAINVAR,COOR,CONNECTb,...
%             TypeElementB,INPUTS_LOC.NameFileMeshLOC_coarse,DATA) ;
    else
        error('Option not implemented')
    end
end
% %else
%     error('This option proved to be inaccurate')
%     [Gb,dR,DOFr,DOFm,AREA,R,DOFA,DOFB,DATAcoarsefine] = ...
%         PERIODIC_BEAMS_2COARSES(a_A,a_B,DOMAINVAR,COOR,CONNECTb,...
%         TypeElementB,INPUTS_LOC.NameFileMeshLOC_coarse,DATA) ;
% end

DOMAINVAR.RigidBodyModes =R ;


