function [WHICH_TRAJ_IS_EACH_SNAP]=  GetTrajectoryEach_Snapshot(trajectoryEACHsnap,nsnaps) ; 
 
if nargin == 0
    load('tmp2.mat')
end

% Identify cluster with training trajectory
% --------------------------------------------
WHICH_TRAJ_IS_EACH_SNAP = zeros(nsnaps,1) ; 
%Ind_Ini_Cluster_Traj = cell(length(TrajectoryEachTimeStep),1) ; 
for itraj = 1:length(trajectoryEACHsnap)
    STEPS_traj = trajectoryEACHsnap{itraj} ; 
     WHICH_TRAJ_IS_EACH_SNAP(STEPS_traj) = itraj ; 
%    [iAAA,iBBB] =   ismember(STEPS_traj,IND_CLUSTERS) ;
%    % Indexes Clusters Trajecetory itraj 
%    IndClustLOC = iBBB(iAAA) ; 
%    IndClustLOC = setdiff(IndClustLOC,1) ; 
%    Ind_Ini_Cluster_Traj{itraj} = [num2str(min(IndClustLOC)),':',num2str(max(IndClustLOC)) ]; 
end