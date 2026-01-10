function  [Ddamp,alphaD,betaD]    = EstimateDdampFromNatFrequencies_EIFEM(OPERFE,DISP_CONDITIONS,DATA) 
%--------------------------------------------------------------------------
%  EstimateDdampFromNatFrequencies_EIFEM
%
%  Estimates the Rayleigh damping coefficients (α, β) such that:
%
%     D_damp = α * M + β * K
%
%  where M is the mass matrix and K is the stiffness matrix. The parameters α and β 
%  are computed to achieve a prescribed damping ratio ξ (e.g., ξ = 0.1) for the first 
%  few natural frequencies of the system.
%
%  This method is applied within the EIFEM framework (Empirical Interscale FEM), adapted
%  for use in both reduced-order (HROM) and full-order models under small-strain dynamics.
%
%  INPUTS:
%    - OPERFE            : structure containing system matrices:
%         > M: global mass matrix
%         > KstiffLINEAR: linear stiffness matrix (assembled)
%    - DISP_CONDITIONS   : contains DOFl, the vector of free DOFs (excluding Dirichlet)
%    - DATA              : optional parameters, including:
%         > CoefficientsDampingAllNatFrequencies: target damping ratio (e.g., 0.05–0.1)
%
%  OUTPUTS:
%    - Ddamp   : damping matrix (Rayleigh type)
%    - alphaD  : coefficient of mass-proportional damping
%    - betaD   : coefficient of stiffness-proportional damping
%
%  METHOD:
%    1. Extract submatrices M_ll and K_ll for free DOFs.
%    2. Solve the generalized eigenvalue problem:
%           Kll * x = λ * M * x   → ω = sqrt(λ)
%    3. Impose target damping ratio for first few frequencies:
%           ξ = 0.5 (α/ω + β*ω)
%       and solve in least-squares sense to determine α and β.
%
%  NOTES:
%    - This method minimizes the damping ratio mismatch in a least-squares sense.
%    - Uses `lsqnonneg` to ensure physical (non-negative) damping coefficients.
%    - Assumes undamped free-vibration modes (no modal damping matrix required).
%
%  ROLE IN EIFEM:
%    - Used when `DATA.EstimateCoefficientsDampingFromNatFrequencies == 1`
%    - Enables dynamic simulations with user-prescribed damping ratio without manually tuning α, β
%
%  AUTHOR:
%    Joaquín A. Hernández Ortega, UPC/CIMNE
%    Campus Nord & Balmes 185, Barcelona
%    Versions: 17–23 May 2024
%    Comments by ChatGPT4, 13-May-2025
%  SEE ALSO:
%    - UndampedFREQ (for solving eigenproblem)
%    - Section 9.3 in dynamic FEM texts (Rayleigh damping)
%    - HROM and NONL examples in 108_EIFEM_metamat directory
%
%--------------------------------------------------------------------------




% Compute Ddamp  =  alpha*M + beta*K, where alpha and beta are estimated
% from a given damping ration   0 <= xiREF <=1
% Adaptation of EstimateDdampFromNatFrequencies.m, developed by 
% JAHO, 17-MAY-2024, Campus Nord, Barcelona, 2024 
% in /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/108_EIFEM_metamat/02_HROMstand.mlx
% Adaptation in /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/108_EIFEM_metamat/04_EIFEM_NONL.mlx
% JAHO, 23-May-2024, Balmes, 185,  Barcelona 
% ----------------------------------------------------
if nargin == 0
    load('tmp3.mat')
end

% Block M_ll 
DOFl = DISP_CONDITIONS.DOFl ; 
M = OPERFE.M(DOFl,DOFl) ; 
% Block K_ll
K = OPERFE.KstiffLINEAR ; % BasisUall'*(OTHER_output.K*BasisUall) ;
Kll = K(DOFl,DOFl) ; 

neig  = size(Kll,1) ; 

 [MODES, FREQ] = UndampedFREQ(M,Kll,neig)  ; 
 
 nfreq_selected = 2; 
 
 FREQ_equat = FREQ(1:nfreq_selected); 
 disp('Natural periods HROM (s)')
 disp(2*pi./FREQ) 
 
 DATA = DefaultField(DATA,'CoefficientsDampingAllNatFrequencies',0.1) ; 
 
      xiDmin = DATA.CoefficientsDampingAllNatFrequencies ;
 
    % Least squares
    b = xiDmin*ones(size(FREQ_equat)); %[xiDmin 1]' ;
    A =  0.5*[1./FREQ_equat FREQ_equat] ; 
    A_all = 0.5*[1./FREQ FREQ] ; 
    
    
    SOL =  lsqnonneg(A,b) ; 
    xi = (A_all*SOL);
    disp(['Damping ratios obained from  imposing the prescribed value for the first ',num2str(nfreq_selected),' frequencis'])
    disp(xi)
    
    betaD = SOL(2) ;
    alphaD = SOL(1) ;
    disp(['betaD = ',num2str(betaD)])
        disp(['alphaD = ',num2str(alphaD)])

%     betaD = 1e-5;
%     disp('borra esto')
    
    Ddamp = alphaD*OPERFE.M + betaD*K ;  
    
    
%     % Damping ratios
%     xiD = 0.5*(alphaD./FREQ + betaD*FREQ) ;
%     aaaa = find(xiD >1);
%     
%     if ~isempty(aaaa)
%         xiD(aaaa) = 1 ;
%     end

 