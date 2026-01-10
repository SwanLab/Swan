function [PhiALL,Phi_non,SNAPdisp_plast_L] =PlasticModesTOLgiven(PhiMaster_lin,SNAPdisp,DOFl,ind_plastic,DATAoffline,Kll_chol,Kll)

if nargin == 0
    load('tmp1.mat')
    DATAoffline.errorDISP = 1e-3 ; 
end

if ~isfield(DATAoffline,'errorDISP')
    error('You shoud specify DATAoffline.errorDISP')
end

TOLsvd = DATAoffline.errorDISP; 

 % We first apply the SVD to the unaltered displacement snapshots
[Phi_non_All_W,SSS,VVV] = SRSVD(SNAPdisp(DOFl,ind_plastic),DATAoffline.errorDISP) ; 
 
 
[Phi_non_All_W,SSS,VVV] = SRSVD(Kll_chol*SNAPdisp(DOFl,ind_plastic),DATAoffline.errorDISP) ; 


% PLASTIC COMPONENT. IT SHOULD BE ORTHOGONAL TO  PhiMaster_lin in the norm
% induced by Kll. The following operator computes the orthogonal projection
[SNAPdisp_plast_L,~,~] =  SprojDEF_operator(PhiMaster_lin,Kll,SNAPdisp(DOFl,ind_plastic)) ;
% Basis matrix for  SNAPdisp_plast_L via weighted SVD. We use randomized
% SVD for efficiency
DATAoffline = DefaultField(DATAoffline,'errorDISP',0) ;
[Phi_non_All_W,SSS,VVV] = SRSVD(Kll_chol*SNAPdisp_plast_L,DATAoffline.errorDISP) ;
% TRUNCATION BASED ON DATAoffline.errorDISP tends to give rather high
% number of plastic modes
% This is why the user should provide an additional input
% DATAoffline.nmodes_PLASTIC
DATAoffline = DefaultField(DATAoffline,'nmodes_PLASTIC',length(SSS)) ;
nmodes_PLASTIC = min(DATAoffline.nmodes_PLASTIC,length(SSS)) ;
Phi_non_All_W = Phi_non_All_W(:,1:nmodes_PLASTIC) ;
Phi_non  = Kll_chol\Phi_non_All_W(:,1:nmodes_PLASTIC) ;
SingValuesPlast = SSS(1:nmodes_PLASTIC)/SSS(1) ;

% Examine total error
PhiALL = [PhiMaster_lin,Phi_non] ;
coeffs = PhiALL'*Kll*SNAPdisp(DOFl,:) ;

ERROR_phi_mat = SNAPdisp(DOFl,:) - PhiALL*coeffs ;
nE2 = sum(sum(ERROR_phi_mat.*(Kll*ERROR_phi_mat))) ;
nD2 = sum(sum(SNAPdisp(DOFl,:).*(Kll*SNAPdisp(DOFl,:)))) ;
totalERROR = sqrt(nE2/nD2) ;
disp(['Total ERROR svd in norm Kll = ',num2str(totalERROR),' for ','nmodes PLASTIC =',num2str(size(Phi_non,2))])


