function CheckAccurayMAT_ECMl2(ECMdata,SNAPfint,W,qLATENT,wORIG)  

if nargin ==0
    load('tmp3.mat')
    close all
end

% Exact integral


w  = feval(ECMdata.wRED.DATA_regress.nameFunctionEvaluate,qLATENT,ECMdata.wRED.DATA_regress) ; 
z = ECMdata.setPoints ; 
bEXACT = cell(1,length(SNAPfint)) ; 
bAPPROX = bEXACT;  
bAPPROX_noreg = bEXACT ; 
% nbE = zeros(size(bEXACT)) ; 
% nbA = nbE; 
% nbANR = nbA ; 
for icluster = 1:length(SNAPfint)
    bEXACT{icluster} = SNAPfint{icluster}'*W ; 
 %   nbE(icluster) = norm(bEXACT{icluster},'fro') ; 
    bAPPROX{icluster} = SNAPfint{icluster}(z,:)'*w(:,icluster) ; 
    bAPPROX_noreg{icluster}   = SNAPfint{icluster}(z,:)'*wORIG(:,icluster) ; 
  %   nbA(icluster) = norm(bAPPROX{icluster},'fro') ; 
   %   nbANR(icluster) = norm(bAPPROX_noreg{icluster},'fro') ; 
end


bEXACT =cell2mat(bEXACT) ; 
bAPPROX = cell2mat(bAPPROX) ; 
bAPPROX_noreg = cell2mat(bAPPROX_noreg) ; 

eee = norm(bEXACT-bAPPROX)/norm(bEXACT) ; 
disp(['Error approx. MAW-ECM, regression (over 1) =',num2str(eee)]) 

eee = norm(bEXACT-bAPPROX_noreg)/norm(bEXACT) ; 
disp(['Error approx. MAW-ECM, no regression (over 1) =',num2str(eee)]) 


qLATENT_show = 1; 
if  qLATENT_show == 0
    qLATENT = 1:size(bEXACT,2) ; 
end

figure(1940)
hold on 
xlabel('qLATENT')
ylabel('Internal forces ')
ncomps = size(SNAPfint{1},2)  ; 
for icomp = 1:ncomps
plot(qLATENT,bEXACT(icomp,:),'DisplayName',['Comp = ',num2str(icomp),'Exact integral, M=',num2str(length(W)),' points'])
plot(qLATENT,bAPPROX(icomp,:),'DisplayName',['Comp = ',num2str(icomp),'Approx. integral, m=',num2str(size(w,1)),' points'])
plot(qLATENT,bAPPROX_noreg(icomp,:),'DisplayName',['Comp = ',num2str(icomp),'Approx. integral no reg.  m=',num2str(size(w,1)),' points'])
end

legend show

