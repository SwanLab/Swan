function [qLATENT_LOC,qL_extended] = EncoderMaxDISP(DATAoffline,SNAP_cluster,DOFl,DATA_evaluateTAU_and_DER,qLATENT_LOC,iloc,MESH)
% ENCODERMAXDISP  Encoder of a scalar latent variable based on the location
%                 of the maximum displacement on a selected surface.
%
%   [qLATENT_LOC,qL_extended] = ENCODERMAXDISP(DATAoffline,SNAP_cluster,DOFl,...
%                       DATA_evaluateTAU_and_DER,qLATENT_LOC,iloc,MESH)
%
%   This routine defines a scalar latent (master) coordinate qMASTER as the
%   x–coordinate of the point where the absolute value of the displacement
%   on a prescribed surface attains its maximum, after a certain snapshot
%   post-processing. This qMASTER is then passed to a user-defined mapping
%   that returns an extended latent vector qL_extended.
%
%   The function is intended to be used as an "encoder" in a manifold-based
%   ROM: for each cluster (or snapshot group), it compresses the displacement
%   snapshots into a single scalar describing the position of the maximum
%   displacement on the chosen boundary surface.
%
%   INPUTS
%   ------
%   DATAoffline :
%       Struct with offline configuration data. Fields used here:
%         .SurfaceMaximumDisplacementLatentVariable
%             Index (integer) of the boundary surface where the maximum
%             displacement is searched. This is used to select the set of
%             nodes from MESH.NODES_FACES.
%         .DirectionMaximumDisplacementLatentVariable
%             Component of the displacement vector to be considered when
%             searching the maximum. Typical values:
%                 1 -> x–direction
%                 2 -> y–direction
%                 3 -> z–direction
%
%   SNAP_cluster :
%       Struct containing snapshot information for the current cluster.
%       The field SNAP_cluster.DISP is expected to be an SVD-like
%       decomposition with fields:
%         .U : matrix of left singular vectors (nodal DOFs x n_sing),
%         .S : vector of singular values (or a diagonal stored as a vector),
%         .V : matrix of right singular vectors (n_sing x n_snapshots).
%       In this encoder, the displacement snapshots are reconstructed (up to
%       a truncation) as
%           SNAPdisp = (U * diag(S)) * V(2:end,:)' ,
%       i.e. using all singular vectors except the first one.
%
%   DOFl :
%       Array with the list of local/global DOFs associated with the cluster.
%       It is included for interface compatibility, but not explicitly used
%       in the current implementation of this encoder.
%
%   DATA_evaluateTAU_and_DER :
%       Struct with information required to build an extended latent vector
%       from the master coordinate. Fields used here:
%         .nameFunctionEvaluate
%             String with the name of a function of the form
%                 [qL_extended,aux] = nameFunctionEvaluate(qMASTER,DATA_evaluateTAU_and_DER)
%             which receives the scalar qMASTER and returns an extended latent
%             vector qL_extended (and optionally additional data).
%
%   qLATENT_LOC :
%       Cell array that stores, for each cluster, the encoded latent variable(s).
%       This function updates the entry qLATENT_LOC{iloc} with the new
%       qMASTER computed for the current cluster.
%
%   iloc :
%       Integer index indicating which entry of qLATENT_LOC corresponds to
%       the current cluster/snapshot group.
%
%   MESH :
%       Mesh data structure. Fields used here:
%         .COOR
%             Matrix of nodal coordinates of size (nnode x ndim).
%         .NODES_FACES
%             Cell array where NODES_FACES{i} contains the node indices of
%             surface/face i. The index
%                 iSURFtop = DATAoffline.SurfaceMaximumDisplacementLatentVariable
%             selects the surface on which the maximum displacement is computed.
%
%   OPTIONAL BEHAVIOUR
%   ------------------
%   If the function is called with no input arguments (nargin == 0), it
%   attempts to load input data from the file 'tmp.mat'. This is meant only
%   for debugging or interactive testing and should not be used in production.
%
%   OUTPUTS
%   -------
%   qLATENT_LOC :
%       Updated cell array of latent coordinates. After the call,
%       qLATENT_LOC{iloc} = qMASTER, where qMASTER contains the x–coordinate
%       (in CoorTOP) of the node at which the absolute displacement on the
%       selected surface and direction reaches its maximum, for each snapshot.
%
%   qL_extended :
%       Extended latent vector returned by the user-supplied evaluation
%       function DATA_evaluateTAU_and_DER.nameFunctionEvaluate. Typically,
%       this collects qMASTER and possibly other derived latent variables
%       required by the ROM.
%
%   ALGORITHM (HIGH-LEVEL)
%   ----------------------
%   1) Identify the top surface (iSURFtop) and displacement direction (idimMAX)
%      from DATAoffline.
%   2) Extract the nodes on that surface (NodesTOP) and build the corresponding
%      DOFs (DofsTOP), as well as their x–coordinates (CoorTOP).
%   3) Reconstruct the displacement snapshots in compressed form using
%          SNAPdisp = (U * diag(S)) * V(2:end,:)' .
%   4) From SNAPdisp, extract the displacements on the selected surface and
%      direction (disTOP).
%   5) For each snapshot, find the node index where |disTOP| is maximal and
%      set qMASTER equal to the corresponding coordinate in CoorTOP.
%   6) Call the user-defined function nameFunctionEvaluate with qMASTER to
%      obtain qL_extended.
%   7) Store qMASTER in qLATENT_LOC{iloc}.
%

if nargin == 0
    load('tmp.mat')
end

iSURFtop = DATAoffline.SurfaceMaximumDisplacementLatentVariable ; 
idimMAX = DATAoffline.DirectionMaximumDisplacementLatentVariable ; 
 NodesTOP = MESH.NODES_FACES{iSURFtop}   ;
   ndim = size(MESH.COOR,2) ;
  DofsTOP = small2large(NodesTOP,ndim) ;
   CoorTOP = MESH.COOR(NodesTOP,1) ;

  
  % coeff =  BasisU(:,1:nREDcoor)'*SNAP_cluster.DISP.U(DOFl,:) ;
%             coeff =bsxfun(@times,coeff',SNAP_cluster.DISP.S)' ;

 SNAPdisp= bsxfun(@times,SNAP_cluster.DISP.U',SNAP_cluster.DISP.S)' ;
 SNAPdisp =  SNAPdisp*SNAP_cluster.DISP.V(2:end,:)' ;
  
    disTOP = SNAPdisp(DofsTOP(idimMAX:ndim:end),:) ; 

  [maxDISP,indXmax] = max(abs(disTOP)) ;
   qMASTER=  CoorTOP(indXmax) ;
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