function  CheckEncoder(BasisU,GmatN,SNAP_cluster,qL_extended,DOFl) 
if nargin == 0
    load('tmp1.mat')
    close all
end


coeffsMasterModes = BasisU'*(GmatN*SNAP_cluster.DISP.U(DOFl,:)) ;  % SVD, U,S,V
coeffsMasterModes =bsxfun(@times,coeffsMasterModes',SNAP_cluster.DISP.S)' ;
qL_FE  = coeffsMasterModes*SNAP_cluster.DISP.V' ;

error_MODES = norm(qL_extended-qL_FE,'fro') ; 
error_MODES = error_MODES/norm(qL_FE,'fro') ;
disp(['error ENCODER = ',num2str(error_MODES)])

figure(1445)
hold on
xlabel('STEP')
ylabel('q')
for i = 1:size(qL_FE,1)
    plot(qL_FE(i,:),'DisplayName',['FE, i=',num2str(i)])
    plot(qL_extended(i,:),'DisplayName',['ENCODER i=',num2str(i)],'LineStyle','--')
end
legend show