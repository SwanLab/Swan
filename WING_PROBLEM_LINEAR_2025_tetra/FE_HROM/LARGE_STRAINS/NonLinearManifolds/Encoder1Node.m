function [qLATENT_LOC,qL_extended] = Encoder1Node(DATAoffline,SNAP_cluster,DOFl,DATA_evaluateTAU_and_DER,qLATENT_LOC,iloc,MESH,DATA_interp)
 
if nargin == 0
    load('tmp.mat')
end


ndim  =size(MESH.COOR,2) ; 

nodeSELECTED = DATA_interp.NodeToFollowAsLatentVariable ; 
dofLOC = DATA_interp.DirectionDisplacementAsLatentVariable ;
DOFlatent = small2large(nodeSELECTED,ndim) ;
DOFlatent = DOFlatent(dofLOC) ; 
% 
% iSURFtop = DATAoffline.SurfaceMaximumDisplacementLatentVariable ; 
% idimMAX = DATAoffline.DirectionMaximumDisplacementLatentVariable ; 
%  NodesTOP = MESH.NODES_FACES{iSURFtop}   ;
%    ndim = size(MESH.COOR,2) ;
%   DofsTOP = small2large(NodesTOP,ndim) ;
%    CoorTOP = MESH.COOR(NodesTOP,1) ;

  
  % coeff =  BasisU(:,1:nREDcoor)'*SNAP_cluster.DISP.U(DOFl,:) ;
%             coeff =bsxfun(@times,coeff',SNAP_cluster.DISP.S)' ;

 SNAPdisp= bsxfun(@times,SNAP_cluster.DISP.U',SNAP_cluster.DISP.S)' ;
 SNAPdisp =  SNAPdisp*SNAP_cluster.DISP.V(2:end,:)' ;
  
%    disTOP = SNAPdisp(DofsTOP(idimMAX:ndim:end),:) ; 

 % [maxDISP,indXmax] = max(abs(disTOP)) ;
   qMASTER=  SNAPdisp(DOFlatent,:) ;
    [qL_extended,~] = feval(DATA_evaluateTAU_and_DER.nameFunctionEvaluate,qMASTER,DATA_evaluateTAU_and_DER) ;
   qLATENT_LOC{iloc} =  qMASTER ;
   
%             
%             coeff =  BasisU(:,1:nREDcoor)'*SNAP_cluster.DISP.U(DOFl,:) ;
%             coeff =bsxfun(@times,coeff',SNAP_cluster.DISP.S)' ;
%             qL  = coeff*SNAP_cluster.DISP.V' ;
%             [qL_extended,~] = feval(DATA_evaluateTAU_and_DER.nameFunctionEvaluate,qL,DATA_evaluateTAU_and_DER) ;
%             qLATENT_LOC{iloc} = qL_extended(idimLAT,:) ;
%             tauNONder = [] ; 
%             
            
%             
%             iSURFtop = 3; 
% DATAoffline = DefaultField(DATAoffline,'SurfaceMaximumDisplacementLatentVariable',iSURFtop) ; 
% idimMAX = 2;  
% DATAoffline = DefaultField(DATAoffline,'DirectionMaximumDisplacementLatentVariable',idimMAX) ; 
% idimMAX = DATAoffline.DirectionMaximumDisplacementLatentVariable ; 
% 
% iSURFtop = DATAoffline.SurfaceMaximumDisplacementLatentVariable ; 
% SNAPdisp = cell2mat(SNAPdisp) ;
% 
% nDOFS = size(SNAPdisp,1) ; 
% NodesTOP = MESH.NODES_FACES{iSURFtop}  ;
% ndim = size(MESH.COOR,2) ;
% DofsTOP = small2large(NodesTOP,ndim) ;
% CoorTOP = MESH.COOR(NodesTOP,1) ;
% disTOP = SNAPdisp(DofsTOP(idimMAX:ndim:end),:) ;
% 
% figure(109)
% hold on
% xlabel('x')
% ylabel('displacement')
% title('Displacement (y) top surface versus x')
% freq = 9 ;
% TIme_select = 1:freq:size(disTOP,2)
% 
% for itime = 1:length(TIme_select)
%     plot(CoorTOP,disTOP(:,TIme_select(itime)),'DisplayName',['FE, istep = ',num2str(TIme_select(itime))]) ;
% end
% 
% [maxDISP,indXmax] = max(abs(disTOP)) ;
% qMASTER=  CoorTOP(indXmax) ; % This play the role of latent variable (master variable)
%  
% figure(190)
% hold on
% xlabel('Time step')
% title('Location Maximum Displacement (latent variable), check it is injective')
% ylabel('q')
% plot(qMASTER)