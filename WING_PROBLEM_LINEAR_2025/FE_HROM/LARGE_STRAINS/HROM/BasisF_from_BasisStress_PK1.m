function SNAPforceS = BasisF_from_BasisStress_PK1(BstRED_l, BasisPone, DATA)
%--------------------------------------------------------------------------
% Computes the reduced internal force snapshots by projecting the stress 
% basis onto the reduced B-operator (BstRED_l). Each block in the result 
% corresponds to the contribution of one displacement mode combined with 
% all stress basis modes. 
%--------------------------------------------------------------------------

% If no input arguments are passed (e.g., for testing), load defaults
if nargin == 0
    load('tmp1.mat')
end

% Ensure the field 'BasisSTRESS_SINGULAR_VALUES' exists in DATA
DATA = DefaultField(DATA, 'BasisSTRESS_SINGULAR_VALUES', []);

% If singular values from SVD are provided, scale the stress basis accordingly
if ~isempty(DATA.BasisSTRESS_SINGULAR_VALUES)
    S = DATA.BasisSTRESS_SINGULAR_VALUES;
    % Normalize with respect to the first singular value
    BasisPone = bsxfun(@times, BasisPone', S / S(1))'; 
end

% Number of displacement modes (columns in BstRED_l)
nDEF = size(BstRED_l, 2); 

% Number of stress basis modes (columns in BasisPone)
nBasisPone = size(BasisPone, 2); 

% Initialize output matrix
SNAPforceS = []; 

% Ensure the flags exist in DATA
DATA = DefaultField(DATA, 'NO_USE_Deformation_gradient_in_Small_Strains', 0); 
DATA = DefaultField(DATA, 'SMALL_STRAIN_KINEMATICS', 0); 

% Determine the number of stress/strain components (nstrainF)
% For small strains: use DATA.MESH.nstrain
% For finite strains: use ndim^2 (e.g., 4 for 2D, 9 for 3D)
if DATA.NO_USE_Deformation_gradient_in_Small_Strains == 1 && DATA.SMALL_STRAIN_KINEMATICS == 1  
    nstrainF = DATA.MESH.nstrain;  
else
    nstrainF = DATA.MESH.ndim^2; 
end

% Loop over all displacement modes
for I = 1:nDEF
    % Initialize the local block of force integrands for displacement mode I
    % Each row corresponds to a Gauss point; each column to a stress basis mode
    SNAPloc = zeros(DATA.MESH.ngausT, nBasisPone); 

    % Loop over all components of the stress tensor (e.g., xx, yy, xy)
    for istrain = 1:nstrainF
        % Extract entries for the current strain component at all Gauss points
        % and multiply with corresponding BstRED_l entries for mode I
        % (This forms a Gauss-wise Hadamard product)
        B_stress = bsxfun(@times, ...
                          BasisPone(istrain:nstrainF:end, :), ...
                          BstRED_l(istrain:nstrainF:end, I)); 
        
        % Accumulate the contribution from this strain component
        SNAPloc = SNAPloc + B_stress;
    end

    % Concatenate the results for this displacement mode
    SNAPforceS = [SNAPforceS, SNAPloc]; 
end
