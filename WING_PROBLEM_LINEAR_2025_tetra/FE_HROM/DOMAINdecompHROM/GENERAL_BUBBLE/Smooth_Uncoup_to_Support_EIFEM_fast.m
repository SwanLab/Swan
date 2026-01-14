function SMOOTH_from_UNCOUP_TO_SUPPORT = Smooth_Uncoup_to_Support_EIFEM_fast(COOR_SUPPORT,IndicesNEW,nnode_UNCOUP)
% 
% To optimize the MATLAB code for creating and manipulating sparse matrices, 
%the primary goal is to avoid populating the sparse matrix within a loop whenever possible. 
%This can be achieved by pre-computing the row indices, column indices, and values that need to be inserted into the matrix,
%and then creating the matrix in one step. Here’s a revised version of your MATLAB code that uses this more efficient approach:

% Number of support nodes
nnode_SUPPORT = size(COOR_SUPPORT, 1);
% Number of uncoupled nodes
 ndim  =size(COOR_SUPPORT,2) ; 
% Extend the index array to ensure it covers all nodes
IndicesNEW_aug = [IndicesNEW; nnode_UNCOUP + 1];

% Prepare to collect indices and values for the sparse matrix
row_indices = [];
col_indices = [];
values = [];

% Generate the row indices, column indices, and values
for inodeSUPP = 1:nnode_SUPPORT
    INDE_LOC = IndicesNEW_aug(inodeSUPP):(IndicesNEW_aug(inodeSUPP + 1) - 1);
    num_elements = length(INDE_LOC);
    
    % Append current indices and values
    row_indices = [row_indices; repmat(inodeSUPP, num_elements, 1)];
    col_indices = [col_indices; INDE_LOC(:)];
    values = [values; repmat(1 / num_elements, num_elements, 1)];
end

% Create the sparse matrix P in one operation
P = sparse(row_indices, col_indices, values, nnode_SUPPORT, nnode_UNCOUP);

% Initialize the larger sparse matrix
SMOOTH_from_UNCOUP_TO_SUPPORT = sparse(ndim*nnode_SUPPORT, ndim*nnode_UNCOUP);

% Assign matrix P to blocks of SMOOTH_from_UNCOUP_TO_SUPPORT
for idim = 1:ndim
    row_range = idim:ndim:(ndim*nnode_SUPPORT);
    col_range = idim:ndim:(ndim*nnode_UNCOUP);
    SMOOTH_from_UNCOUP_TO_SUPPORT(row_range, col_range) = P;
end
% ```
% 
% ### Key Changes Made
% 1. **Pre-computation of Indices and Values**: The indices for rows and columns along with their corresponding values are computed outside of the sparse matrix creation. This method is more memory-efficient and significantly faster than modifying a sparse matrix iteratively within a loop.
% 
% 2. **Bulk Insertion into Sparse Matrix**: By using the `sparse` function with precomputed row indices, column indices, and values, the matrix `P` is created in a single step. This technique avoids the overhead of repeated memory allocation associated with incrementally adding elements to a sparse matrix.
% 
% 3. **Efficient Block Assignment**: In the second loop, instead of operating element-by-element, entire blocks of the matrix are filled by assigning the smaller matrix `P` directly to slices of the larger matrix, leveraging MATLAB’s efficient handling of such operations on sparse matrices.
% 
% This approach should provide a substantial improvement in performance, particularly for large matrices, by minimizing the overhead associated with the dynamic structure of sparse matrix storage in MATLAB.
% 
% ndim  = size(COOR_SUPPORT,2) ; 
% nnode_SUPPORT = size(COOR_SUPPORT,1) ;
%  P = sparse(nnode_SUPPORT,nnode_UNCOUP) ;
% 
% IndicesNEW_aug = [IndicesNEW;nnode_UNCOUP] ;
% for inodeSUPP = 1:nnode_SUPPORT
%     INDE_LOC = IndicesNEW_aug(inodeSUPP):(IndicesNEW_aug(inodeSUPP+1)-1) ;
%     P(inodeSUPP,INDE_LOC) =  1/length(INDE_LOC)  ;
% end
% SMOOTH_from_UNCOUP_TO_SUPPORT = sparse(ndim*nnode_SUPPORT,ndim*nnode_UNCOUP) ;
% for idim = 1:ndim
%     SMOOTH_from_UNCOUP_TO_SUPPORT(idim:ndim:end,idim:ndim:end) = P ;
% end