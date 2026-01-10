function [rDEF ]= AmplitudeReactions(DATAROM,MESH1D,a,fextBEAMr,ndim) ;

if nargin ==0
    load('tmp1.mat')
end

% if length(DATAROM) == 1
%     % Only slices of one single type
%     [rDEF ]= AmplitudeReactions_onlybeam(DATAROM{1},MESH1D,a,fextBEAMr,ndim)  ;
%     
% else
    
    % Slices and joints
    [rDEF ]= AmplitudeReactions_jbeam(DATAROM,MESH1D,a,fextBEAMr,ndim)  ;
    
% end