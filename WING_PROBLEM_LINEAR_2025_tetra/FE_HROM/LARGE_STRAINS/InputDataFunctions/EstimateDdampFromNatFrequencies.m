function  [Ddamp,alphaD,betaD]    = EstimateDdampFromNatFrequencies(OPERHROM,OTHER_output,BasisUall,DISP_CONDITIONS,DATAHROM) 
% Compute Ddamp  =  alpha*M + beta*K, where alpha and beta are estimated
% from a given damping ration   0 <= xiREF <=1
% JAHO, 17-MAY-2024, Campus Nord, Barcelona, 2024 
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/108_EIFEM_metamat/02_HROMstand.mlx
% ----------------------------------------------------
if nargin == 0
    load('tmp3.mat')
end

% Block M_ll 
DOFl = DISP_CONDITIONS.DOFl ; 
M = OPERHROM.M(DOFl,DOFl) ; 
% Block K_ll
K = BasisUall'*(OTHER_output.K*BasisUall) ;
Kll = K(DOFl,DOFl) ; 

neig  = size(Kll,1) ; 

 [MODES, FREQ] = UndampedFREQ(M,Kll,neig)  ; 
 
 nfreq_selected = 2; 
 
 FREQ_equat = FREQ(1:nfreq_selected); 
 disp('Natural periods HROM (s)')
 disp(2*pi./FREQ) 
 
 DATAHROM = DefaultField(DATAHROM,'CoefficientsDampingAllNatFrequencies',0.1) ; 
 
      xiDmin = DATAHROM.CoefficientsDampingAllNatFrequencies ;
 
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
    
%     betaD = 1e-5;
%     disp('borra esto')
    
    Ddamp = alphaD*OPERHROM.M + betaD*K ;  
    
    
%     % Damping ratios
%     xiD = 0.5*(alphaD./FREQ + betaD*FREQ) ;
%     aaaa = find(xiD >1);
%     
%     if ~isempty(aaaa)
%         xiD(aaaa) = 1 ;
%     end

 