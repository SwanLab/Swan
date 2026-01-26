function [MESH,CONVERSION_FROM_OLD_TO_NEW] = RemoveInteriorPointHexa27(MESH)
if nargin == 0
    load('tmp.mat')
end

% THE GOAL HERE IS TO CONSTRUCT A COLUMN MATRIX "CONVERTION_NEW_OLD", IN
% WHICH YOU INPUT A LIST OF NODES IN THE OLD MESH  (WHICH INCLUDES THE INTERIOR NODES)
% AND IT RETURNS THE CORRESPONDING LABELLING IN THE NEW 

%  MESH.CN
%     11     1     9    22    23    10    21    27     3     2    12    16    17     4    13    25    14    15    24    26     5     7     6    18    20    19  8   8
NodesRemove = MESH.CN(:,end) ;
% NodesRemove =
% 
%      8
CN = MESH.CN(:,1:end-1) ;
% For instance
%CN =

%    11     1     9    22    23    10    21    27     3     2    12    16    17     4    13    25    14    15    24    26     5     7     6    18    20    19
NodesRemain = setdiff(1:size(MESH.COOR,1),NodesRemove) ;
% CONVERSION INDEX VECTOR 
CONVERSION_FROM_OLD_TO_NEW  = zeros(size(MESH.COOR,1),1) ; 
CONVERSION_FROM_OLD_TO_NEW(NodesRemain) = 1:length(NodesRemain) ; 

%%% 


%MESH.COOR_INTPOINT = MESH.COOR(NodesRemove,:); 
MESH.COOR = MESH.COOR(NodesRemain,:) ;
% NodesRemain =
% 
%      1     2     3     4     5     6     7     9    10    11    12    13    14    15    16    17    18    19    20    21    22    23    24    25    26    27

%MESH.CN = RenumberConnectivities(CN,1:length(NodesRemain)) ;
CNprev =CN(:) ; 
CNnew = CONVERSION_FROM_OLD_TO_NEW(CNprev) ; 
MESH.CN = reshape(CNnew,[],26) ; 




% MESH.CN
% 
% ans =
% 
%     10     1     8    21    22     9    20    26     3     2    11    15    16     4    12    24    13    14    23    25     5     7     6    17    19    18
% MESH.CNb = RenumberConnectivities(MESH.CNb,1:length(NodesRemain)) ;
CNprev =MESH.CNb(:) ; 
CNnew = CONVERSION_FROM_OLD_TO_NEW(CNprev) ; 
MESH.CNb = reshape(CNnew,[],9) ; 
% % Before 
% MESH.CNb
% 
% ans =
% 
%     11     1     9    22     3     2    12    16     5
%      1    10    21     9     4    15    13     2     6
%     11    23    10     1    17    14     4     3     7
%      9    21    27    22    13    24    25    12    18
%     23    27    21    10    26    24    15    14    19
%     11    22    27    23    16    25    26    17    20
% AFTER
% MESH.CNb
% 
% ans =
% 
%     10     1     8    21     3     2    11    15     5
%      1     9    20     8     4    14    12     2     6
%     10    22     9     1    16    13     4     3     7
%      8    20    26    21    12    23    24    11    17
%     22    26    20     9    25    23    14    13    18
%     10    21    26    22    15    24    25    16    19




if  isfield(MESH,'NODES_FACES') &&  ~isempty(MESH.NODES_FACES{1})
    for i = 1:length(MESH.NODES_FACES)
        MESH.NODES_FACES{i} = CONVERSION_FROM_OLD_TO_NEW(MESH.NODES_FACES{i}) ;
        
        % Before 
%         MESH.NODES_FACES{i}'
% 
% ans =
% 
%      9    12    13    18    21    22    24    25    27
  % After 
%   MESH.NODES_FACES{i}'
% 
% ans =
% 
%      8    11    12    17    20    21    23    24    26
  
        
    end
end



if isfield(MESH,'NODES_LINES') && ~isempty(MESH.NODES_LINES)
    for i = 1:length(MESH.NODES_LINES)
        MESH.NODES_LINES{i} = CONVERSION_FROM_OLD_TO_NEW(MESH.NODES_LINES{i}) ;
    end
end

if isfield(MESH,'NODES_DOM') &&  ~isempty(MESH.NODES_DOM)
    for i = 1:length(MESH.NODES_DOM)
        MESH.NODES_DOM{i} = CONVERSION_FROM_OLD_TO_NEW(MESH.NODES_DOM{i}) ;
    end
end