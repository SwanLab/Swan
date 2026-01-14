function [CtangM,  stressNP1, EPnp1,  alphaNP1,  sigmayNP1 ] = ...
    SmallStrainJ2Plasticity_LOCAL(eNP1,PMAT,EPn,sigmayN,alphaN,DATA)
%--------------------------------------------------------------------------
% function [CtangM, stressNP1, EPnp1, alphaNP1, sigmayNP1] = ...
%     SmallStrainJ2Plasticity_LOCAL(eNP1, PMAT, EPn, sigmayN, alphaN, DATA)
%
% PURPOSE:
%   Implements the update of stresses and internal variables for an
%   isotropic J2 small-strain elastoplastic constitutive model with linear
%   isotropic hardening, following the return-mapping algorithm of 
%   Simo & Hughes (see *Computational Inelasticity*, 1998).
%
%   The function is fully vectorized over Gauss points, making it suitable
%   for finite element implementations where stresses and strains are stored
%   in a global (stacked) format.
%
% ALGORITHM:
%   1. Compute elastic trial strains from total strain minus plastic strain.
%   2. Form trial deviatoric and spherical stresses.
%   3. Evaluate the J2 yield condition at each Gauss point.
%   4. Split Gauss points into elastic and plastic sets:
%        - Elastic: accept trial stress, tangent = elastic tensor.
%        - Plastic: perform radial return mapping:
%             * Compute plastic multiplier Δγ.
%             * Update stresses, plastic strains, isotropic hardening 
%               variable (alpha), and yield stress.
%             * Update consistent tangent operator.
%
% INPUTS:
%   eNP1     - [ngausT × 1] vector of total strains at time step n+1
%   PMAT     - structure with material parameters, containing:
%                 • mu        : shear modulus at each Gauss point
%                 • kappa     : bulk modulus
%                 • Hmodulus  : isotropic hardening modulus
%                 • Celas     : elastic constitutive matrix (sparse format)
%                 • identVOL  : volumetric projection tensor
%                 • identV    : deviatoric projection tensor
%   EPn      - [ngausT × 1] plastic strain vector at step n
%   sigmayN  - [ngaus × 1] yield stress at step n
%   alphaN   - [ngaus × 1] accumulated plastic strain at step n
%   DATA     - structure with solver options and mesh info, containing:
%                 • MESH.nstrain : number of strain components per Gauss pt
%                 • MESH.ngaus   : number of Gauss points
%                 • CALC_CTANG   : flag (1/0) to compute consistent tangent
%                 • kiter        : current Newton iteration
%                 • FirstIterationAtWhichYieldingIsActivated : 
%                      iteration threshold to activate plasticity check
%
% OUTPUTS:
%   CtangM     - [ngausT × nstrain] consistent tangent operator matrix
%   stressNP1  - [ngausT × 1] updated Cauchy stress components at step n+1
%   EPnp1      - [ngausT × 1] updated plastic strain vector at step n+1
%   alphaNP1   - [ngaus × 1] updated accumulated plastic strain
%   sigmayNP1  - [ngaus × 1] updated yield stress (isotropic hardening)
%
% NOTES:
%   - The function applies a *radial return algorithm* in stress space.
%   - Vectorized indexing is used through the COMP struct for efficiency.
%   - For elastic points, stresses are directly computed from the trial
%     state and tangent = Celas.
%   - For plastic points, stresses and tangent are corrected consistently
%     with the return mapping scheme.
%   - The check DATA.kiter > FirstIterationAtWhichYieldingIsActivated is
%     used to delay plasticity activation to avoid inconsistencies in the
%     first Newton iteration.
%
% REFERENCE:
%   J.C. Simo & T.J.R. Hughes, *Computational Inelasticity*, Springer, 1998.
% JAHO, comments updated by ChatGPT, 16-08-2025
%--------------------------------------------------------------------------



% This is a copy of /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/FE_HROM/FE_CODE/NONLINEAR_Dyn/UpdateStresses_J2_FE2D.m
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/FE_HROM/FE_CODE/NONLINEAR_Dyn/UpdateStresses_J2_FE2D.m

if nargin == 0
    load('tmp.mat')  
end

DATA = DefaultField(DATA,'FirstIterationAtWhichYieldingIsActivated',1) ;%  See definition 


ngaus = size(alphaN,1); ngausT = size(eNP1,1); nstrain = DATA.MESH.nstrain; %   ngausT/ngaus ;
ngausE =DATA.MESH.ngaus ;
COMP.c1 = 1:nstrain:ngausT ;  %  Indices of Component x of stresses and strains
COMP.c2 = 2:nstrain:ngausT ;  % Indices of Component y of stresses and strains
COMP.c3 = 3:nstrain:ngausT ; % Indices of Component xy of stresses and strains
COMP.c4 = 4:nstrain:ngausT ;  % Indices of Component z of stresses and strains
% muA = PMAT.muA ; kappaA = PMAT.kappaA ; kappa = PMAT.kappa ;
% mu = PMAT.mu ;
% Initializations
EPnp1 = EPn ;
alphaNP1 = alphaN ;
sigmayNP1 = sigmayN;
%CtangNP1_B = sparse(size(BBfNW,1),size(BBfNW,2));
%CtangNP1_B= BBfNW ;
CtangM = PMAT.Celas ;
%stressNP1 = zeros(size(EPnp1));
% 1. Trial step
% -------------
% Trial elastic strains
eTRelas = eNP1-EPn ;
% Deviatoric part and spherical parts
tr_elasTR = TracePlaneStrain(eTRelas,COMP) ;  % Trace
dev_elasTR = eTRelas -1/3*tr_elasTR ;
% Deviatoric trial stress
muA =reshape([PMAT.mu'; PMAT.mu'; PMAT.mu'; PMAT.mu'],ngausT,1) ;
sTR = StrainMultiplStress(dev_elasTR,2*muA,COMP) ;
% Spheric trial stress
kappaA =reshape([PMAT.kappa'; PMAT.kappa'; PMAT.kappa'; PMAT.kappa'],ngausT,1) ;
pTR = kappaA.*tr_elasTR ;

% 2. Check yield condition
%-------------------------
normSTR = sqrt(sTR(COMP.c1).^2 +sTR(COMP.c2).^2 +2*sTR(COMP.c3).^2    +sTR(COMP.c4).^2 ) ;
fTR = normSTR -sqrt(2/3)*sigmayN ;
indPLAS = find(fTR>0) ;  % Plastic gauss points
indELAS = 1:ngaus ; indELAS(indPLAS) = [] ;  % Elastic gauss points
%----------------------------------------------
% --- POINTS BEHAVING ELASTICALLY
%----------------------------------------------
%indELASm = INDEXcomponents(indELAS,nstrain)  ;
% By default
%stressNP1(indELASm)  = sTR(indELASm) + pTR(indELASm) ;
stressNP1 = sTR + pTR ;

%  Product CtangST*BBnwCALC_CTANG
% ----------------------
% METHOD_CTANG = 2 ;
% if  CALC_CTANG == 1
%     if METHOD_CTANG ==1
%         CtangNP1_B(indELASm,:)  = CtangB_elas(PMAT,BBfNW,indELASm,nstrain,indELAS) ;
%     end
% end
%----------------------------------------------
% --- POINTS BEHAVING PLASTICALLY
%----------------------------------------------
%dbstop('64')
niter_elast = DATA.FirstIterationAtWhichYieldingIsActivated ; 
if ~isempty(indPLAS) && DATA.kiter>niter_elast   
    % This second condition  (DATA.kiter>niter_elast, where niter_elast =1) is used so that the return mapping algorithm is not activated in the very first iteration
    % JAHO, 30-Oct-2023. The above condition may lead to inconsistencies in the iterative algorithm. More specifically, if
    % it happens that the elastic trial state is in global equilibrium, but
    % violates the yield condition. This has been reported in 
    % /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/104_EIFEM_plast2D/08_IndividualTRAINbub.mlx
    % This is why we have introduce DATA.FirstIterationAtWhichYieldingIsActivated
    
    
    indPLASm = INDEXcomponents(indPLAS,nstrain)  ;
    % Plastic multiplier
    pmult = fTR(indPLAS)./(2*PMAT.mu(indPLAS) + 2/3*PMAT.Hmodulus(indPLAS)) ;
    ngausTg =length(indPLASm) ;
    COMPg.c1 = 1:nstrain:ngausTg ;  %  Indices of Component x of stresses and strains
    COMPg.c2 = 2:nstrain:ngausTg ;  % Indices of Component y of stresses and strains
    COMPg.c3 = 3:nstrain:ngausTg ; % Indices of Component xy of stresses and strains
    COMPg.c4 = 4:nstrain:ngausTg ;  % Indices of Component z of stresses and strains
    pmultA =   reshape([pmult'; pmult'; pmult'; pmult'],ngausTg,1) ; % At All components
    % Internal variables
    alphaNP1(indPLAS) = alphaN(indPLAS) + sqrt(2/3)*pmult;
    sigmayNP1(indPLAS) = sigmayN(indPLAS) + PMAT.Hmodulus(indPLAS).*(alphaNP1(indPLAS)- alphaN(indPLAS));
    %
    % Normal along the plastic flow direction
    normSTRg =  normSTR(indPLAS) ;
    normSTRgA =   reshape([normSTRg'; normSTRg'; normSTRg'; normSTRg'],ngausTg,1) ;  % All components
    normalNP = sTR(indPLASm)./normSTRgA ;
    % Deviatoric stress
    sNP = sTR(indPLASm) - 2*muA(indPLASm).*pmultA.*normalNP ;
    % Total stress
    stressNP1(indPLASm)  =   pTR(indPLASm)  + sNP ;
    % Plastic strains
    devdEP =   StressMultiplStrain(normalNP,pmultA,COMPg) ;
    EPnp1(indPLASm) =  EPn(indPLASm) + devdEP ;
    % Product C*BBfNW
    if  DATA.CALC_CTANG == 1
        %      if METHOD_CTANG ==1
        %         CtangNP1_B(indPLASm,:)  = CtangB_plas(PMAT,BBfNW,indPLASm,nstrain,indPLAS,normSTRg,pmult,normalNP) ;
        %    else
        HmodulusG = PMAT.Hmodulus(indPLAS);
        HmodulusG =reshape([HmodulusG'; HmodulusG'; HmodulusG'; HmodulusG'],length(indPLASm),1) ;
        CtangPLAS  =  Ctang_plasMethod2(nstrain,normSTRgA,pmultA...
            ,normalNP,muA(indPLASm),kappaA(indPLASm),HmodulusG,PMAT.identVOL(indPLASm,:),PMAT.identV(indPLASm,:)) ;
        CtangM(indPLASm,:) = CtangPLAS ;
        %   end
        
    end
    
% elseif  ~isempty(indPLAS) && iter<=niter_elast
%     stressNP1(indPLAS)  = sTR(indPLAS) + pTR(indPLAS) ;
%     %   if METHOD_CTANG ==1
%     %     CtangNP1_B(indPLAS,:)  = CtangB_elas(PMAT,BBfNW,indPLASm,nstrain,indPLAS) ;
%     %  end
end

% if  CALC_CTANG ==1
%     %dbstop('106')
%     if ~isempty(WEIGHTS)
%         for istrain = 1:size(CtangM,2)
%             CtangM(:,istrain) =  CtangM(:,istrain).*WEIGHTS ;
%         end
%         
%     end
%     
%     
%     Ctang =     ConvertCmatSparseMatrix(CtangM,PMAT.indI,PMAT.indJ) ;
%   %  CtangNP1_B = Ctang*BBfNW ;
% end
% % 
% 
% TEST_SOLUTION = 1 ; 
% if TEST_SOLUTION ==1 
%     stress = Ctang*eNP1 ;
% end

