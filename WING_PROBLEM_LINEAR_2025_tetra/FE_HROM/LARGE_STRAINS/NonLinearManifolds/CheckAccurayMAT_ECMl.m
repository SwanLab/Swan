function CheckAccurayMAT_ECMl(ECMdata,SNAPfint,W,qLATENT,wALL)  

if nargin ==0
    load('tmp.mat')
end

% Exact integral


w  = feval(ECMdata.wRED.DATA_regress.nameFunctionEvaluate,qLATENT,ECMdata.wRED.DATA_regress) ; 
z = ECMdata.setPoints ; 
bEXACT = cell(length(SNAPfint),1) ; 
bAPPROX = bEXACT;  
bAPPROX_noreg = bEXACT;  

nbE = zeros(size(bEXACT)) ; 
nbA = nbE; 
for icluster = 1:length(SNAPfint)
    bEXACT{icluster} = SNAPfint{icluster}'*W ; 
    nbE(icluster) = norm(bEXACT{icluster},'fro') ; 
    bAPPROX{icluster} = SNAPfint{icluster}(z,:)'*w(:,icluster) ; 
        bAPPROX_noreg{icluster} = SNAPfint{icluster}(z,:)'*wALL(:,icluster) ; 

    
     nbA(icluster) = norm(bAPPROX{icluster},'fro') ; 
end


bEXACT =cell2mat(bEXACT) ; 
bAPPROX = cell2mat(bAPPROX) ; 
eee = norm(bEXACT-bAPPROX)/norm(bEXACT) ; 
disp(['Error approx. MAW-ECM (over 1) =',num2str(eee)]) 

figure(194)
hold on 
xlabel('qLATENT')
ylabel('Integral ')
plot(qLATENT,nbE,'DisplayName',['Exact integral, M=',num2str(length(W)),' points'])
plot(qLATENT,nbA,'DisplayName',['Approx. integral, m=',num2str(size(w,1)),' points'])
legend show

