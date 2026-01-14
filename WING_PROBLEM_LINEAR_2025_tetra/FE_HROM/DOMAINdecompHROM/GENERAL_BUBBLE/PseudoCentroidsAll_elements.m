function PSEUDO_centroidsELEM = PseudoCentroidsAll_elements(CN,COORref)
% --------------------------------------------------------------------------
% PseudoCentroidsAll_elements
%
% This function computes the pseudo-centroids of all elements in a mesh.
% The pseudo-centroid of an element is defined as the arithmetic mean of 
% the coordinates of its nodes.
%
% INPUTS:
%   CN       : (nelem x nnode_per_elem) connectivity matrix. Each row 
%              contains the global node indices for one element.
%   COORref  : (nnode x ndim) matrix of nodal coordinates. Each row is the 
%              coordinate vector of a node in the reference configuration.
%
% OUTPUT:
%   PSEUDO_centroidsELEM : (nelem x ndim) matrix. Each row is the pseudo-
%                           centroid of the corresponding element.
%
% OPTIONS:
%   VECTORIZED = 1 : Use vectorized implementation for improved performance.
%              = 0 : Use classical loop-based implementation for clarity.
%
% JAHO/ChatGPT-4, May 5th 2025, Campus Nord UPC, Barcelona
%--------------------------------------------------------------------------

VECTORIZED = 1;
if VECTORIZED == 0
    nnodesPERELEMENT = size(CN,2) ;
    PSEUDO_centroidsELEM = zeros(size(CN,1),size(COORref,2)) ;
    for ielemREF = 1:size(CN,1)
        % Nodal coordinates
        nodesLOC = CN(ielemREF,:) ;
        COORoneELEM = COORref(nodesLOC,:) ;
        % Pseudo centroid
        PScentroidELEM = sum(COORoneELEM,1)/nnodesPERELEMENT ;
        PSEUDO_centroidsELEM(ielemREF,:) = PScentroidELEM ;
    end
else
    % VECTORIZED VERSION (ChatGPT 4)
    nnodesPERELEMENT = size(CN, 2);
    % Gather all nodal coordinates for each element
    COORallELEM = reshape(COORref(CN',:), nnodesPERELEMENT, [], size(COORref,2));
    COORallELEM = permute(COORallELEM, [2 1 3]); % Shape: (nelements, nnodesPERELEMENT, ndim)
    % Compute mean (pseudo-centroid) along the node dimension
    PSEUDO_centroidsELEM = mean(COORallELEM, 2);
    PSEUDO_centroidsELEM = squeeze(PSEUDO_centroidsELEM); % Final shape: (nelements, ndim)
    
    
end