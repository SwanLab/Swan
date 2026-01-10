function INTERPOLATED_BASES =  InterpolationSeveralSubspaces...
    (ModesMatrix,FACTOR_SCALES_ALL,FACTOR_TEST,SCALE_FACTOR,DATAINPUT)

if nargin == 0
    load('tmp.mat')
end

% Page 11 of An Interpolation Method for Adapting Reduced-Order
% -------------------------------------------------------------

% STEP 0:
% A point Si0 of the manifold is chosen as a reference and origin point for the interpolation.
iref = 1 ;

% STEP 1
% The tangent space TSi0 and those points among {Si}NR−1
%i=0 which lie in a sufficiently small
% neigborhood of Si0 are now considered. More specifically, each point Si that is sufficiently close to Si0
% is mapped to a matrix 􀀀i representing a point i of TSi0 using the logarithm map LogSi0 . This can
% be written as


PHI_0 =  ModesMatrix{iref} ;

[UU,SS,VV] =  SVDT(PHI_0) ;
TOLLOC = 1e-4 ;
if abs(sum(SS)-length(SS))> TOLLOC
    APPLY_SVD_ON_MODES = 1 ;
    PHI_0 = UU  ;
else
    APPLY_SVD_ON_MODES = 0 ; % modes are already orthogonal....
end

GAMMA = {} ;
for   ipoint = 1:length(ModesMatrix)
    if  ipoint == iref
        GAMMA{ipoint} = zeros(size( ModesMatrix{iref} )) ;  % REference point in the tangent space... 
    else
        PHI_1 =  ModesMatrix{ipoint} ;
        if APPLY_SVD_ON_MODES == 1
            [PHI_1,SS,VV] =  SVDT(PHI_1) ;
        end
        % Scalar product between the two basis matrices
        SP = PHI_0'*PHI_1;
        % Projection of PHIp onto the orthogonal complement of PHI
        OC = PHI_1- PHI_0*(PHI_0'*PHI_1) ;
        % Matrix over which the SVD is to be applied
        X = OC*inv(SP) ;
        % SVD of X
        [U,S,V] = SVDT(X) ;
        % arc-Tangent of S
        ANGLE_0 = atan(S) ;
        
        U = bsxfun(@times,U',ANGLE_0)' ;
        GAMMA{ipoint} = U*V' ;
    end
end


% Step 2. Each entry of the matrix 􀀀NR associated with the target operating point NR is computed
% by interpolating the corresponding entries of the matrices {􀀀i} 2 RNf×N associated with the operating
% points {i}. The choice of interpolation method depends on the number of physical parameters
% contained in each operating point, Np. When Np = 1, a univariate — typically, a Lagrange-type — interpolation
% method is chosen. Otherwise, a multivariate interpolation scheme (see for example [42,43])
% is chosen.

GAMMA_INTERP = InterpolationMatricesDirect(FACTOR_SCALES_ALL,GAMMA,FACTOR_TEST,SCALE_FACTOR,DATAINPUT) ;

% Step 3. The matrix 􀀀NR representing NR 2 TSi0 is mapped to a subspace SNR on the Grassmann
%manifold spanned by a matrix 	NR using the exponential map ExpSi0 . This can be written as

%[U,S,V] = SVDT(GAMMA_INTERP) ;
[U,S,V] = svd(GAMMA_INTERP,0) ;
S = diag(S);

% Term proportional to PHI_0
cosTERM = cos(S) ;
TERM_0 = PHI_0*V ;
TERM_0 = bsxfun(@times,TERM_0',cosTERM)' ;
% Term proportional to PHI_1 (its  orthogonal complement to PHI_0)
sinTERM = sin(S) ;
TERM_1 = bsxfun(@times,U',sinTERM)' ;

INTERPOLATED_BASES = TERM_0 + TERM_1 ;
