function [Hhred,CmacroROM,CmacroHROM,MSG ]= ReducedHomogOperatorFLUC(Z,Wred,Bred,BstD,nstrain,Hred,Cglo_W,wSTs,BasisS)


% Hyperreduced B matrix
% --------------------------
Zstr = small2large(Z,nstrain);  % Set of indexes associated to the selected Gauss Points
Bhred = Bred(Zstr,:) ;  % B-reduced matrix evaluated at the selected Gauss points
% B(z,:)*D 
BdRED = BstD(Zstr,:) ;
% Elasticity matrix at points Z (multiplied by the FE weights)
Cred_W = Cglo_W(Zstr,Zstr) ;   % For implementation reasons, this matrix contains the FE integration weightws
%
% Product Cred_W*Bhred
Cred_Bhred_W  = Cred_W*Bhred ;  %
Cred_BdRED_W  = Cred_W*BdRED ;  %

% Product B*Wred -Initialization
Bhred_Wred = zeros(size(Bhred)) ;
% PRoduct Cred*Bhred - Initialization
Cred_Bhred = zeros(size(Cred_Bhred_W)) ;
Cred_BdRED = zeros(size(Cred_BdRED_W)) ;
for istrain = 1:nstrain
    INDS = istrain:nstrain:length(Zstr) ;
    Bhred_Wred(INDS,:) = bsxfun(@times,Bhred(INDS,:),Wred) ;
    Cred_Bhred(INDS,:) =     bsxfun(@times,Cred_Bhred_W(INDS,:),1./wSTs(Z)) ;
    Cred_BdRED(INDS,:) =     bsxfun(@times,Cred_BdRED_W(INDS,:),1./wSTs(Z)) ;
end
% Hyperreduced matrix
Khred = Bhred_Wred'*Cred_Bhred ;%
fhred = Bhred_Wred'*Cred_BdRED ; 
Hhred = Khred\fhred ;

errorKF = norm(Hred-Hhred,'fro')/norm(Hred,'fro')*100 ; 
MSG = {} ; 
MSG{end+1} = ['Error in approximating matrix Hred=   inv(Kred)*fRED (%) = ',num2str(errorKF)] ; 
disp(MSG{end} )


% In terms of elasticity matrix (macro)
% --------------------------------------
% Reduced-order model   
%  
q = - Hred ; 
strain =-Bred*Hred  + BstD  ; 
stress = Cglo_W*strain ; 
CmacroROM = zeros(nstrain, nstrain  ) ; 
VOL = sum(wSTs) ; 
 for istrain = 1:nstrain
     COMP = istrain:nstrain:size(Cglo_W,1) ;
     CmacroROM(istrain,:) = ones(size(wSTs'))*stress(COMP,:)/VOL ;
 end
 
 
 % HyperReduced-order model   
%  -------------------------------------------
% Reconstruction matrix 
coeff = (BasisS(Zstr,:)'*BasisS(Zstr,:))\BasisS(Zstr,:)' ;
ReconsStresses = BasisS*coeff ; 
% Homogeneized reconstruction stresses 
HOMOGENIZATION_OPERATOR = zeros(nstrain,size(ReconsStresses,2)) ;
%
 for istrain = 1:nstrain
     COMP = istrain:nstrain:size(ReconsStresses,1) ;
     HOMOGENIZATION_OPERATOR(istrain,:) = wSTs'*ReconsStresses(COMP,:)/VOL ;
 end
strain = -Bred(Zstr,:)*Hhred + BstD(Zstr,:) ; 
stress = Cglo_W(Zstr,Zstr)*strain ; 
 for istrain = 1:nstrain
     COMP = istrain:nstrain:size(stress,1) ;
     stress(COMP,:) = bsxfun(@times,stress(COMP,:),1./wSTs(Z)) ;
 end

CmacroHROM =HOMOGENIZATION_OPERATOR*stress ; 


errorKF = norm(CmacroROM-CmacroHROM,'fro')/norm(CmacroROM,'fro')*100 ; 
MSG{end+1} = ['Error in approximating Cmacro  (%) = ',num2str(errorKF)] ; 
disp(MSG{end} )
 
