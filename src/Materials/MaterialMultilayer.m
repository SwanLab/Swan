classdef MaterialMultilayer < handle
    
    properties (Access = public)
        ndim
        E      % Young's moduli per layer [nLayers x 3] or [nLayers x 1]
        nu     % Poisson's ratios per layer [nLayers x 3] or [nLayers x 1]
        G      % Shear moduli per layer [nLayers x 3] or [nLayers x 1]
        rotation  % Ply angles per layer [nLayers x 1] (degrees)
        materialType  % 'ISOTROPIC' or 'ORTHOTROPIC'
        z_interfaces  % Layer interface z-coordinates
        h      % Layer thicknesses
        nLayers % Number of layers
    end

    methods (Access = public)
        
        function obj = MaterialMultilayer(cParams)
            obj.init(cParams);
        end

        function C = evaluate(obj, xV)
            % evaluate - Evaluate constitutive tensor at integration points
            %
            % Input:
            %   xV - Coordinates of evaluation points
            %
            % Output:
            %   C - Constitutive tensor (3 x 3 x 3 x 3 x nGauss x nElem)
            
            nGauss = size(xV, 2);
            nElem  = size(xV, 3);
            
            
            % mid-plane
            z_eval = zeros(nGauss, nElem);
            
            C = zeros(3, 3, 3, 3, nGauss, nElem);
            
            for e = 1:nElem
                for g = 1:nGauss
                    kLayer = obj.findLayer(z_eval(g, e));
                    % Get 6x6 Voigt matrix for this layer
                    C_matrix = obj.createConstitutiveMatrixForLayer(kLayer);
                    
                    % Apply rotation if orthotropic
                    if strcmp(obj.materialType, 'ORTHOTROPIC') && obj.rotation(kLayer) ~= 0
                        theta = deg2rad(obj.rotation(kLayer));
                        C_matrix = obj.rotateConstitutiveMatrix(C_matrix, theta);
                    end
                    % Convert Voigt matrix to 4th order tensor
                    C_tensor = obj.voigtToTensor(C_matrix);                    
                    C(:,:,:,:,g,e) = C_tensor;
                end
            end
        end
        
        function C_tensor = getConstitutiveTensorForLayer(obj, kLayer, rotated)
            % getConstitutiveTensorForLayer - Get 4th order tensor for specific layer
            %
            % Input:
            %   kLayer  - Layer number
            %   rotated - (optional) Apply rotation if true (default: true)
            %
            % Output:
            %   C_tensor - Constitutive tensor (3x3x3x3)
            
            if nargin < 3
                rotated = true;
            end
            
            % Get 6x6 Voigt matrix
            C_matrix = obj.createConstitutiveMatrixForLayer(kLayer);
            % Apply rotation if requested
            if rotated && strcmp(obj.materialType, 'ORTHOTROPIC') && obj.rotation(kLayer) ~= 0
                theta = deg2rad(obj.rotation(kLayer));
                C_matrix = obj.rotateConstitutiveMatrix(C_matrix, theta);
            end
            % Convert to 4th order tensor
            C_tensor = obj.voigtToTensor(C_matrix);
        end
        
        function C_ps = planeStressReduction(obj,C)
            % planeStressReduction - Reduce a 3D constitutive tensor to plane stress
            %
            % For shells/plates, the condition sigma_33 = 0 is enforced.
            % This is applied via Schur-complement condensation on the z-direction (index 3):
            %
            %   C*_abcd = C_abcd - C_ab33 * C_3333^{-1} * C_33cd
            %
            % where a,b,c,d in {1,2}.
            % The out-of-plane shear terms (i3),(3j) are NOT affected and are kept as-is.

            C_ps = C;
            C33 = C(3,3,3,3);
            if abs(C33) < 1e-30
                return  % Nothing to condense (already degenerate)
            end
            for i = 1:3
                for j = 1:3
                    for k = 1:3
                        for l = 1:3
                            % Only condense in-plane normal components
                            if (i ~= 3 && j ~= 3 && k ~= 3 && l ~= 3)
                                C_ps(i,j,k,l) = C(i,j,k,l) - C(i,j,3,3)*C(3,3,k,l)/C33;
                            end
                        end
                    end
                end
            end
        end
        
    end
    
    methods (Access = public)
        
        function init(obj, cParams)
            obj.ndim = cParams.ndim;
            obj.materialType = cParams.materialType;            
            if isa(cParams.E, 'ConstantFunction')
                obj.E = cParams.E.constant;
            else
                obj.E = cParams.E;
            end
            
            if isa(cParams.nu, 'ConstantFunction')
                obj.nu = cParams.nu.constant;
            else
                obj.nu = cParams.nu;
            end
            
            if isfield(cParams, 'h')
                if isa(cParams.h, 'ConstantFunction')
                    obj.h = cParams.h.constant;
                else
                    obj.h = cParams.h;
                end
            else
                obj.h = 1; 
            end
            obj.nLayers = length(obj.h);           
            if strcmp(obj.materialType, 'ORTHOTROPIC')
                if isa(cParams.G, 'ConstantFunction')
                    obj.G = cParams.G.constant;
                else
                    obj.G = cParams.G;
                end
                
                if isfield(cParams, 'rotation')
                    obj.rotation = cParams.rotation;
                else
                    obj.rotation = zeros(obj.nLayers, 1);
                end
            else
                % Isotropic: compute G
                obj.G = obj.E ./ (2 .* (1 + obj.nu));
                obj.rotation = zeros(obj.nLayers, 1);
            end            
            H_total = sum(obj.h);
            obj.z_interfaces = cumsum([0; obj.h(:)]) - H_total/2;
        end
        
        function kLayer = findLayer(obj, z)
            % Find which layer z belongs to
            kLayer = find(z >= obj.z_interfaces(1:end-1) & ...
                          z <= obj.z_interfaces(2:end), 1);
            if isempty(kLayer)
                % Out of bounds: use closest layer
                if z < obj.z_interfaces(1)
                    kLayer = 1;
                else
                    kLayer = obj.nLayers;
                end
            end
        end
        
        function C_matrix = createConstitutiveMatrixForLayer(obj, kLayer)
            % Create 6x6 Voigt matrix for a specific layer
            % Voigt order: [11, 22, 33, 23, 13, 12]
            if strcmp(obj.materialType, 'ISOTROPIC')
                if size(obj.E, 2) == 1
                    E_k = obj.E(kLayer);
                    nu_k = obj.nu(kLayer);
                else
                    E_k = obj.E(kLayer, 1);
                    nu_k = obj.nu(kLayer, 1);
                end
                
                lambda = E_k * nu_k / ((1 + nu_k) * (1 - 2*nu_k));                
                C_matrix = lambda * [1-nu_k, nu_k,   nu_k,   0,             0,             0;
                                     nu_k,   1-nu_k, nu_k,   0,             0,             0;
                                     nu_k,   nu_k,   1-nu_k, 0,             0,             0;
                                     0,      0,      0,      (1-2*nu_k),    0,             0;
                                     0,      0,      0,      0,             (1-2*nu_k),    0;
                                     0,      0,      0,      0,             0,             (1-2*nu_k)];
                
            else
                E1 = obj.E(kLayer, 1);
                E2 = obj.E(kLayer, 2);
                E3 = obj.E(kLayer, 3);
                
                nu12 = obj.nu(kLayer, 1);
                nu13 = obj.nu(kLayer, 2);
                nu23 = obj.nu(kLayer, 3);
                
                G12 = obj.G(kLayer, 1);  % G12 (1-2 plane)
                G13 = obj.G(kLayer, 2);  % G13 (1-3 plane)
                G23 = obj.G(kLayer, 3);  % G23 (2-3 plane)
                
                S = zeros(6, 6);
                S(1,1) = 1/E1;  S(2,2) = 1/E2;  S(3,3) = 1/E3;
                S(1,2) = -nu12/E1; S(2,1) = S(1,2);
                S(1,3) = -nu13/E1; S(3,1) = S(1,3);
                S(2,3) = -nu23/E2; S(3,2) = S(2,3);
                S(4,4) = 1/G23;  % 23 shear
                S(5,5) = 1/G13;  % 13 shear
                S(6,6) = 1/G12;  % 12 shear
                
                C_temp = inv(S);          
                C_matrix = C_temp;
                % C_matrix(4:6, 4:6) = C_temp(4:6, 4:6) * 2;
            end
        end
        
        function C_tensor = voigtToTensor(obj, C_matrix)
            % voigtToTensor - Convert 6x6 Voigt matrix to 3x3x3x3 tensor
            %
            % Voigt order: [11, 22, 33, 23, 13, 12] (STANDARD)
            %
            % Input:
            %   C_matrix - 6x6 Voigt matrix 
            %
            % Output:
            %   C_tensor - 3x3x3x3 constitutive tensor
            
            % Initialize tensor
            C_tensor = zeros(3, 3, 3, 3);
            
            % Voigt index mapping (STANDARD): [11, 22, 33, 23, 13, 12]
            voigt_map = [1,1;  % 1 → (1,1)
                         2,2;  % 2 → (2,2)
                         3,3;  % 3 → (3,3)
                         2,3;  % 4 → (2,3)
                         1,3;  % 5 → (1,3)
                         1,2]; % 6 → (1,2)
            
            for p = 1:6
                i = voigt_map(p, 1);
                j = voigt_map(p, 2);
                
                for q = 1:6
                    k = voigt_map(q, 1);
                    l = voigt_map(q, 2);
                    
                    % Get value from Voigt matrix
                    C_val = C_matrix(p, q);                   
                    C_tensor(i,j,k,l) = C_val;
                    C_tensor(j,i,k,l) = C_val;  % Symmetry in (i,j)
                    C_tensor(i,j,l,k) = C_val;  % Symmetry in (k,l)
                    C_tensor(j,i,l,k) = C_val;  % Both symmetries
                    C_tensor(k,l,i,j) = C_val;  % Major symmetry
                    C_tensor(l,k,i,j) = C_val;
                    C_tensor(k,l,j,i) = C_val;
                    C_tensor(l,k,j,i) = C_val;
                end
            end
        end
        
        function C_matrix = tensorToVoigt(obj, C_tensor)
            % tensorToVoigt - Convert 3x3x3x3 tensor to 6x6 Voigt matrix
            % Voigt order: [11, 22, 33, 23, 13, 12] (STANDARD)
            
            voigt_map = [1,1; 2,2; 3,3; 2,3; 1,3; 1,2];
            C_matrix = zeros(6, 6);
            
            for p = 1:6
                i = voigt_map(p, 1);
                j = voigt_map(p, 2);
                
                for q = 1:6
                    k = voigt_map(q, 1);
                    l = voigt_map(q, 2);
                    
                    C_val = C_tensor(i, j, k, l);
                    
                    % Apply Voigt factor (multiply by 2 for shear)
                    if p > 3, C_val = C_val * 2; end
                    if q > 3, C_val = C_val * 2; end
                    
                    C_matrix(p, q) = C_val;
                end
            end
        end
        
        function C_rotated = rotateConstitutiveMatrix(obj, C, theta)
            % Rotate Voigt matrix using transformation matrix
            T = obj.getTransformationMatrix(theta);
            C_rotated = T * C * T.';
        end
        
        function T = getTransformationMatrix(obj, theta)
            % Transformation matrix for STANDARD Voigt notation
            % Voigt order: [11, 22, 33, 23, 13, 12]
            
            c = cos(theta);
            s = sin(theta);
            c2 = c^2;
            s2 = s^2;
            s2theta = sin(2*theta);
            sc = s*c;
            
            % Standard transformation matrix
            T = [c2,   s2,  0,  0,          0,         -s2theta;
                 s2,   c2,  0,  0,          0,          s2theta;
                 0,    0,   1,  0,          0,          0;
                 0,    0,   0,  c,          s,          0;
                 0,    0,   0, -s,          c,          0;
                 sc,  -sc,  0,  0,          0,          c2-s2];
        end

        
    end
end