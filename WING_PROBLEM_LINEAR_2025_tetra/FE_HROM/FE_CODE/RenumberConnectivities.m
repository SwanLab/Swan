function CNrenumbered = RenumberConnectivities(CN,NODES)
%--------------------------------------------------------------------------
% FUNCTION: RenumberConnectivities
%
% PURPOSE:
%   Renumbers a connectivity matrix `CN` based on a new node ordering `NODES`.
%   It returns a new connectivity matrix `CNrenumbered` in which each global node
%   number in `CN` is replaced with a local index according to its position in `NODES`.
%
% USAGE:
%   CNrenumbered = RenumberConnectivities(CN, NODES)
%
% INPUT:
%   - CN     : (nElem x nNodePerElem) connectivity matrix with global node IDs.
%   - NODES  : (nNodes x 1) vector containing the new local node IDs,
%              in the order corresponding to the local numbering.
%              Typically, NODES = unique(CN(:)) for a subdomain.
%
% OUTPUT:
%   - CNrenumbered : (nElem x nNodePerElem) connectivity matrix with local
%                    node IDs based on `NODES`.
%
% ALGORITHM:
%   - Vectorized implementation:
%     1. Flattens the connectivity matrix into a column.
%     2. Sorts the vector and records index mappings (forward and backward).
%     3. Locates the positions of unique node IDs using `unique`.
%     4. Assigns new IDs from `NODES` based on the sorted unique values.
%     5. Inverts the sorting to restore original element order and reshapes.
%
%   - Fallback (non-vectorized) version:
%     Loops through each unique node in the input, finds occurrences in CN,
%     and assigns a new index from 1 to nNodes. Less efficient but easier to debug.
%
% IMPORTANT NOTES:
%   - The function assumes that `length(NODES) == length(unique(CN(:)))`.
%   - If errors arise, check if the mesh (especially boundary surfaces)
%     has redundant or mismatched node labels.
%
% EXAMPLE:
%   CN = [90 80; 80 60];        % Global IDs
%   NODES = [60 80 90];         % New ordering
%   CNrenumbered = RenumberConnectivities(CN, NODES);
%   % CNrenumbered would then be [3 2; 2 1]
%
% AUTHOR:
%   Joaquín A. Hernández Ortega
%   UPC-CIMNE
%   Version: 11-Apr-2024
%   Comments by ChatGPT4, 13-May-2025
%--------------------------------------------------------------------------

if nargin == 0
    %     % -----------------------------------------------
    %     % Inputs: Connectivity matrix with old numbering
    %     %------------------------------------------------
    %     CN = [90 80
    %         80 60 ] ;
    %     % New numbering so that.
    %     NODES = [1 2 3 4] ;
    %     % length(NODES) must be equal to length(unique(CN(:)))
    %     % Output: Connectivity matrix done according to the new numbering
    %     %
    %     VECTORIZED_VERSION = 0 ;
    
    load('tmp2.mat')
    VECTORIZED_VERSION = 0;
end

VECTORIZED_VERSION =1;

if  VECTORIZED_VERSION == 0
    CNrenumbered =   RenumemberCONNECT_serial(CN,NODES) ;
else
    %%%%
    % This is like having a list of numbers
    CN_col = CN(:) ;
    % Sorting
    [CNsort,Iforward,Iback] = sortJAHO(CN_col) ;
    
    [CNunique LastUnique] = unique(CNsort,'last');
    [CNfirst FirstUnique] = unique(CNsort,'first');
    
    CN_ren = zeros(size(CN_col));
    % JAHO, 11-Apr-2024
    % If the progam fails at this point, check that the surfaces of the
    % geometric model has been correctly labelled (no redundancies...)
    
    for i=1:length(CNunique)
        CN_ren(FirstUnique(i):LastUnique(i)) =NODES(i);
    end
    
    % Now we undo the operations
    % --------------------------
    CNrenumbered = CN_ren(Iback) ;
    CNrenumbered = reshape(CNrenumbered,[],size(CN,2));
    
end

end

function CNnew =   RenumemberCONNECT_serial(CN,NODES) ;
disp('Non-vectorized renumering')
ListNodes = unique(CN(:)) ;
CNnew = zeros(size(CN(:))) ;
for inode = 1:length(ListNodes)
    inodeOLD = ListNodes(inode) ;
    inodeNEW = inode  ;
    IND =  find(CN(:)==inodeOLD) ;
    CNnew(IND) = inodeNEW ;
end
ncol = size(CN,2) ;
CNnew = reshape(CNnew,[],ncol) ;

end

