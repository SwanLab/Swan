function AssemblyPeriodCond(nods,nodm)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Assembly of vector of prescribed displacements and G
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %------------------
nDOFr = 2*length(nods) + length(nods_w) ;  % Number of slave DOFs 
nDOFm = 2*length(nodm) + length(nodm_w) ;  % Number of master DOFm

% -----------
% SLAVES DOFs 
% -----------
[DOFr dR ROWSslv] = AssemDOFslv_PERIOD_TB(nDOFr,nods,nods_w,uB,uB_w) ; 
% -----------
% MASTER DOFs 
% -----------
[DOFm ROWSmst] = AssemDOFmst_PERIOD_TB(nDOFm,nodm,nodm_w) ; 
% Matrix G 
% --------
G= sparse(nDOFr,nDOFm) ;
%%% u, direction 
G(ROWSslv.u,ROWSmst.u) = Gi ; 
%%% v, direction 
G(ROWSslv.v,ROWSmst.v) = Gi ; 
%%% w, direction 
G(ROWSslv.w,ROWSmst.w) = Gi_w ; 

%%%% Now we sort DOFr and DOFm in ascending order
[DOFr iSLV] = sort(DOFr,'ascend') ; 
[DOFm iMST] = sort(DOFm,'ascend') ; 
% Therefore 
dR(iSLV) = dR ; 
G(iSLV,iMST) =G ; 

 


