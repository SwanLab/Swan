function [BasisU,FREQ] = ModesFrequenciesSmallStrains(OPERFE,DOFl,DATA,MESH,MATPRO,nMODES,dINI,SNAPdisp,SNAPvel)

if nargin == 0
    load('tmp2.mat')
end

disp('Computing BasisU as natural modes corresponding to K and M')
ndofALL = size(OPERFE.M,1) ;
dLOC = zeros(ndofALL,1) ;

% 2. Deformation gradient at all Gauss points
FgradST = OPERFE.Bst*dLOC + repmat(OPERFE.IDENTITY_F,1,size(dLOC,2)) ;
% 3. Green-Lagrante strains at all Gauss points
GLSTRAINS = StrainGreenLagrange(FgradST,DATA.MESH.ndim) ;
% 4. 2nd Piola-Kirchhoff stresses at all Gauss Points

[PK2STRESS,celasST ]= PK2stress_Constitutive_Model(GLSTRAINS,MATPRO,DATA,FgradST) ;


K = KstiffLargeStrains(OPERFE,PK2STRESS,FgradST,DATA.MESH.ndim,celasST) ;



SORT_AMPLI = 0 ;

if SORT_AMPLI == 0
    [BasisUprueba,FREQ] = UndampedFREQ(OPERFE.M(DOFl,DOFl),K(DOFl,DOFl),nMODES)  ;

else
    
    
    nMODES_all = 3*nMODES ; 

[BasisUprueba,FREQ] = UndampedFREQ(OPERFE.M(DOFl,DOFl),K(DOFl,DOFl),nMODES_all)  ;
    
    AMPLITUDES= BasisUprueba'*OPERFE.M(DOFl,DOFl)*dINI(DOFl) ;
    
    AMPLITUDES= AMPLITUDES/max(abs(AMPLITUDES))*100 ;
    
    [AMPLITUDES,indSORTED] = sort(abs(AMPLITUDES),'descend') ;
    
    FREQ = FREQ(indSORTED(1:nMODES)) ;
    BasisUprueba = BasisUprueba(:,indSORTED(1:nMODES)) ;
end
disp('Natural frequencies')
disp(num2str(FREQ'))

[BasisU] = SVDT(BasisUprueba) ; 

Mred = BasisU'*OPERFE.M(DOFl,DOFl)*BasisU ;
Kred = BasisU'*K(DOFl,DOFl)*BasisU ;

 [BasisUprueba_red,FREQ_red] = UndampedFREQ(Mred,Kred,nMODES)  ;
 
 
 SNAPdisp = cell2mat(SNAPdisp) ; 
 
 ERROR_SNAPdisp = SNAPdisp-BasisU*(BasisU'*SNAPdisp) ; 
 
 ERROR_snapN = norm(ERROR_SNAPdisp,'fro')/norm(SNAPdisp,'fro') ; 
 
 disp(['Error approximation SNAPdisp with natural modes ( over 1)=',num2str(ERROR_snapN)])
 
  SNAPvel = cell2mat(SNAPvel) ; 
 
 ERROR_SNAPvel = SNAPvel-BasisU*(BasisU'*SNAPvel) ; 
 
 ERROR_snapN = norm(ERROR_SNAPvel,'fro')/norm(SNAPvel,'fro') ; 
 
 disp(['Error approximation SNAPvel with natural modes ( over 1)=',num2str(ERROR_snapN)])
 


%  Sdisp = [] ; Vdisp = [] ;