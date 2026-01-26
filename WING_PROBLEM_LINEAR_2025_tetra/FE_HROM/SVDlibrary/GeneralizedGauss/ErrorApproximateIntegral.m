function IntApproxSVD = ErrorApproximateIntegral(A,Lambda,INTexac,W,DATA_ECM,zDECM,wDECM,DATA)

if nargin == 0
    load('tmp.mat')
end

if ~iscell(A)
    Vw = Lambda'*A ;
    IntApproxSVD = (Lambda*Vw)'*sqrt(W) ; 
    INTdecm = A(zDECM,:)'*(wDECM./sqrt(W(zDECM))) ; % High-fidelity integral of all snapshots
else
    IntApproxSVD = cell(size(A)) ;
      INTdecm = cell(size(A)) ;
    disp('Computing   integrals projected matrices ')
    for imat = 1:length(A)
        disp(['iblock = ',num2str(imat)])
        if isnumeric(A{imat})
            Vw = Lambda'*A{imat} ;
            IntApproxSVD{imat} =  (Lambda*Vw)'*sqrt(W) ;
            INTdecm{imat} =  A{imat}(zDECM,:)'*(wDECM./sqrt(W(zDECM))) ;
        else
            disp(['Retrieving from memory ...'])
            SSS = load(A{imat}) ;
            disp([' ... Done'])
            fff = fieldnames(SSS) ;
            %   Ai = SSS.(fff{1}) ; SSS = [] ;
            Vw = Lambda'*SSS.(fff{1}) ;
            IntApproxSVD{imat} = (Lambda*Vw)'*sqrt(W) ;
             INTdecm{imat} =  SSS.(fff{1})(zDECM,:)'*(wDECM./sqrt(W(zDECM))) ;
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

