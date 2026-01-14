function CheckAccurayMAT_ECM(ECMdata,SNAPfint,OPERFE,qLATENT)  

if nargin ==0
    load('tmp.mat')
end

% Exact integral


w  = feval(ECMdata.wRED.DATA_regress.nameFunctionEvaluate,qLATENT,ECMdata.wRED.DATA_regress) ; 
z = ECMdata.setPoints ; 
bEXACT = cell(length(SNAPfint),1) ; 
bAPPROX = bEXACT;  
for icluster = 1:length(SNAPfint)
    bEXACT{icluster} = SNAPfint{icluster}'*OPERFE.wSTs ; 
    bAPPROX{icluster} = SNAPfint{icluster}(z,:)'*w(:,icluster) ; 
end


bEXACT =cell2mat(bEXACT) ; 
bAPPROX = cell2mat(bAPPROX) ; 
eee = norm(bEXACT-bAPPROX)/norm(bEXACT) ; 
disp(['Error approx. MAW-ECM (over 1) =',num2str(eee)]) 
