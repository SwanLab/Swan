function PERMUT = PermutationConnectivities(TypeElement,nnodeE)
%--------------------------------------------------------------------------
%  PermutationConnectivities
%
%  Returns a list of valid permutations of element connectivities for a
%  given element type (`TypeElement`) and number of nodes (`nnodeE`), such
%  that the geometric Jacobian determinant remains positive in each case.
%
%  This is crucial when assembling finite element matrices using pre-trained
%  operators that are defined in a reference parent domain. The mapping must
%  preserve orientation and avoid negative Jacobians, especially under
%  transformations and rotations (see large deformation handling).
%
%  The permutations are generated according to:
%   - Line elements (2-node)
%   - Quadrilaterals (4-node and 8-node)
%   - Hexahedra (8-node, 20-node, and 26-node)
%   - User-defined cases (no permutations applied)
%
%  For hexahedral elements, dedicated permutation functions are called:
%   - `PermutationConnHexahedra8`
%   - `PermutationConnHexahedra20`
%   - `PermutationConnHexahedra26`
%
%  These functions are expected to return all consistent permutations that
%  preserve positive Jacobian orientation under affine mappings.
%
%  REFERENCES:
%    - Detailed usage and examples in:
%      /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/
%        101_MULTI2D_2023/06_3D_27_HEXAHEDRA/01_8nodeHEX.mlx
%
%  INPUTS:
%    - TypeElement : string ('Linear', 'Quadrilateral', 'Hexahedra', etc.)
%    - nnodeE      : number of nodes per element (e.g., 2, 4, 8, 20, 26)
%
%  OUTPUT:
%    - PERMUT : cell array containing permutations of node indices
%
%  Author:
%    Joaquín A. Hernández Ortega, UPC/CIMNE
%    Balmes 185, Barcelona
%    Version: 25-March-2023
%    Comments by ChatGPT4, 13-May-2025
%
%--------------------------------------------------------------------------




% Given TypeElement and nnodeE (number of nodes per element), this function returns the possible permutations
% of the connectivity of the element so that the Jacobian of the
% transformation is always positive
% See further information in       
% /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/101_MULTI2D_2023/06_3D_27_HEXAHEDRA/01_8nodeHEX.mlx

% JAHO, 25-March-2023
if nargin == 0
    load('tmp.mat')
end

    % POSSIBLE PERMUTATIONS CONNECTIVITIES
    switch TypeElement
        case 'Linear'
            if nnodeE ==2
                PERMUT = {[1,2],[2,1]} ; 
            else
                error('Option not implemented')
            end
        case 'Quadrilateral'
            if nnodeE ==4
                PERMUT = {[1:4],[2:4,1],[3:4,1:2],[4,1:3]};
                % disp('Temporal amendment ! Change it!')
                %PERMUT = {[1:4]} ;
            elseif  nnodeE==8
                PERMUT = {[1:4,5:8],[2:4,1,6:8,5],[3:4,1:2,7:8,5:6],[4,1:3,8,5:7]};
            else
                error('Option not implemented')
            end
        case 'Hexahedra'
            if nnodeE ==8
                % See 
                % /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/101_MULTI2D_2023/06_3D_27_HEXAHEDRA/01_8nodeHEX.mlx
                 PERMUT = PermutationConnHexahedra8 ; 
                 
            elseif  nnodeE ==20
                PERMUT = PermutationConnHexahedra20 ; 
                
            elseif nnodeE == 26 
                PERMUT = PermutationConnHexahedra26 ; 
            else
                 error('Option not implemented')
              
            end
        case 'GIVEN_BY_USER'     
            PERMUT = {[1:nnodeE]} ; 
            
        otherwise
            error('Option not implemented')
    end
    