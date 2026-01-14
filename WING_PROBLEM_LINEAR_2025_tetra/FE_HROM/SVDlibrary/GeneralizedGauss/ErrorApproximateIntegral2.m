function IntApproxSVD = ErrorApproximateIntegral2(A,Lambda,W,DATA_ECM,zDECM,wDECM,DATA)

if nargin == 0
    load('tmp.mat')
end
INTexac = DATA.ExactIntegral ; 
if ~iscell(A)
    Vw = Lambda'*(bsxfun(@times,A,W)) ;
    IntApproxSVD = (Lambda*Vw)'*W ; 
    INTdecm = A(zDECM,:)'*wDECM ;  % High-fidelity integral of all snapshots
else
    IntApproxSVD = cell(size(A)) ;
      INTdecm = cell(size(A)) ;
    disp('Computing   integrals projected matrices ')
    for imat = 1:length(A)
        disp(['iblock = ',num2str(imat)])
        if isnumeric(A{imat})
            Vw = Lambda'*(bsxfun(@times,A{imat},W)) ; 
            IntApproxSVD{imat} =  (Lambda*Vw)'*W ;
            INTdecm{imat} =  A{imat}(zDECM,:)'*wDECM ;
        else
            disp(['Retrieving from memory ...'])
            SSS = load(A{imat}) ;
            disp([' ... Done'])
            fff = fieldnames(SSS) ;
            %   Ai = SSS.(fff{1}) ; SSS = [] ;
            Vw = Lambda'*(bsxfun(@times,SSS.(fff{1}),W)) ;
            IntApproxSVD{imat} = (Lambda*Vw)'*W ;
             INTdecm{imat} =  SSS.(fff{1})(zDECM,:)'*wDECM ;
            %   Ai = bsxfun(@times,Ai,sqrt(W)) ;
            %   disp(['Savaing again in  memory ... (multiplied by sqrt(W))'])
            %  save(A{imat},'Ai') ;
            %    disp([' ... Done'])
        end
    end
    IntApproxSVD = cell2mat(IntApproxSVD') ;
     INTdecm = cell2mat(INTdecm') ;
end

ErrorINtAPPROXsvd = norm(INTexac-IntApproxSVD)/norm(INTexac)*100 ;
disp('----------------------------------------------------------------------------------------')
disp(['Integration Error associated to the SVD =',num2str(ErrorINtAPPROXsvd),' % (for a SVD tolerance =',num2str(DATA_ECM.TOL_SVD_A)  ,' )'])
disp('----------------------------------------------------------------------------------------')

 IntegrationError = INTexac - INTdecm ;
 IntegrationError = norm(IntegrationError)/norm(INTexac)*100;
 disp(['Actual integration error using DECM= ',num2str(IntegrationError),' % (prescribed tolerance fint =',num2str(DATA_ECM.TOL_SVD_A*100), '%'])
% 


% if ~iscell(A)
%     INTdecm = A(zDECM,:)'*(wDECM./sqrt(W(zDECM))) ; % High-fidelity integral of all snapshots
% else
%     INTdecm = cell(size(A)) ;
%     for imat = 1:length(A)
%         INTdecm{imat} =A{imat}(zDECM,:)'*wDECM ;
%     end
%     INTdecm = cell2mat(INTdecm') ;
% end
% IntegrationError = INTexac - INTdecm ;
% IntegrationError = norm(IntegrationError)/norm(INTexac)*100;
% disp(['Actual integration error using DECM= ',num2str(IntegrationError),' % (prescribed tolerance fint =',num2str(DATA_ECM.TOL_SVD_A*100), '%'])
% 

