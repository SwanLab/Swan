function [PhiALL,Phi_non,SNAPdisp_plast_L] =...
    PlasticModesNMODESgiven2(PhiMaster_lin,SNAPdisp,DOFl,ind_plastic,DATAoffline,Kll_chol,Kll)
if nargin == 0
    load('tmp1.mat')
    DATAoffline.nmodes_PLASTIC = 50 ; 
    
end

[Phi_non_All_W,SSS,VVV] = SRSVD(Kll_chol*SNAPdisp(DOFl,ind_plastic),DATAoffline.errorDISP) ;
DATAoffline = DefaultField(DATAoffline,'nmodes_PLASTIC',length(SSS)) ;
nmodes_PLASTIC = min(DATAoffline.nmodes_PLASTIC,length(SSS)) ;
 Phi_non  = Kll_chol\Phi_non_All_W(:,1:nmodes_PLASTIC) ;


[Phi_non_orth,~,~] =  SprojDEF_operator(PhiMaster_lin,Kll,Phi_non) ;

 DATALOCa.RELATIVE_SVD = 1 ; 
  DATALOCa.Mchol = Kll_chol ; 
  DATALOCa.TOL = 1e-10 ; 

[Phi_non_orth, S,~] = WSVDT(Phi_non_orth, [], DATALOCa)  ; 
 
 
% Examine total error
PhiALL = [PhiMaster_lin,Phi_non_orth] ;
coeffs = PhiALL'*Kll*SNAPdisp(DOFl,:) ;

ERROR_phi_mat = SNAPdisp(DOFl,:) - PhiALL*coeffs ;
nE2 = sum(sum(ERROR_phi_mat.*(Kll*ERROR_phi_mat))) ;
nD2 = sum(sum(SNAPdisp(DOFl,:).*(Kll*SNAPdisp(DOFl,:)))) ;
totalERROR = sqrt(nE2/nD2) ;
disp(['Total ERROR svd in norm Kll = ',num2str(totalERROR),' for ','nmodes PLASTIC =',num2str(size(Phi_non,2))])


