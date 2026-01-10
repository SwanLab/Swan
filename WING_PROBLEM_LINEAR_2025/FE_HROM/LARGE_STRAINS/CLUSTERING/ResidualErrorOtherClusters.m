function  ErrorResidualOtherClusters = ResidualErrorOtherClusters(iCL,TransClustDATA,VAR,DATAoffline,DATA_cl,...
    OPERFE_cl,MATPRO_cl,DISP_CONDITIONS_cl,istep,Fbody_cl,Ftrac_cl)

if nargin == 0
    load('tmp1.mat')
end

disp('***********************************************')
disp(['Evaluating residual remaining clusters '])
disp('***********************************************')


q = VAR.DISP(DISP_CONDITIONS_cl{iCL}.DOFl) ;  % Current reduced coordinates
nq2  = norm(q)^2 ;  % Norm current reduced coordinates


ErrorResidualOtherClusters = zeros(size(OPERFE_cl)) ; 
 ErrorTransitionALL =  ErrorResidualOtherClusters ; % = sqrt((nq2 - norm(qNEW)^2)/nq2) ;

for jcluster = 1:length(OPERFE_cl)
  
    VAR.DISP = zeros(size(VAR.DISP)) ; 
    DOFl = DISP_CONDITIONS_cl{jcluster}.DOFl ; 
    %if jcluster~=iCL
    VAR = BoundaryConditionsImpose(VAR,DISP_CONDITIONS_cl,Fbody_cl,Ftrac_cl,jcluster,istep) ; 
    if jcluster == iCL
        qNEW = q ; 
    else
    qNEW = TransClustDATA.TransMatrix{jcluster,iCL}*q ;
    end
    VAR.DISP(DOFl) = qNEW ;     
    [~,~,~,VAR.RESID,Fint,~,~] = ResidualFromDisplacements(OPERFE_cl{jcluster},VAR.DISP,MATPRO_cl{jcluster},DATA_cl{jcluster},VAR.FEXT) ; 
     [ErrorResidualOtherClusters(jcluster)]=  CheckConvergenceLSTR_nomessage(Fint(DOFl),VAR.FEXT(DOFl),VAR.RESID(DOFl),0,1e20) ;
     ErrorTransitionALL(jcluster) = sqrt((nq2 - norm(qNEW)^2)/nq2) ; 
       disp(['k = ',num2str(jcluster),'   norm(res)  =',num2str(ErrorResidualOtherClusters(jcluster)),'  error trans = ',num2str(ErrorTransitionALL(jcluster))])
    %end
  
end
