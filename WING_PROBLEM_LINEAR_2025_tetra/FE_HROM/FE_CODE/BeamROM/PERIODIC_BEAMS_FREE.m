% function [G,uBAR,DOFr,DOFm,AREA,R,DOFA,DOFB] = PERIODIC_BEAMS_FREE(a_A,a_B,DOMAINVAR,COOR,CONNECTb,TypeElementB,...
%     NameFileMeshCoarse,DATA)
% 
% if nargin == 0
%     load('tmp2.mat')
% end
% % 
% % DATA = DefaultField(DATA,'FLUCTUATION_MODES_BENDING',[]) ;
% % Ub = DATA.FLUCTUATION_MODES_BENDING ;
% 
% ndim = 3;
% da = a_A-a_B ;
% 
% 
% %-------------------------------------------------------
% % FINER MESH
% %%%% FACE 1
% iface=1 ;
% nodesfA = DOMAINVAR.NODES_faces12{1,iface} ;
% if  max(nodesfA-sort(nodesfA)) ~=0
%     error(['nodesfA should be sorted in ascending order'])
% end
% DOFA = small2large(nodesfA,ndim) ;   % DOFS face 1
% % Geometric mass matrix, centroid
% COOR_FACE = COOR(nodesfA,:) ; % Coordinates of this face
% 
% [CentroidFA,AREA,Mst] =CentroidGeometricMassMatrixNEW(COOR,nodesfA,CONNECTb{1},TypeElementB) ;
% 
% %[CentroidFB,AREA,Mst] =CentroidGeometricMassMatrixNEW(DATA_REFMESH.COOR,nodesfB,DATA_REFMESH.CONNECTb{iface},TypeElementB) ;
% 
% %
% COORrelA = bsxfun(@minus,COOR_FACE',CentroidFA')'; % Coordinates relative to centroid
% R = ConstructBasisRigidBody(COORrelA) ;
% 
% % FACE 2
% iface=2 ;
% nodesfB = DOMAINVAR.NODES_faces12{end,iface} ;   % They are already paired,
% DOFB = small2large(nodesfB,ndim) ; % Set of DOFs face B
% % -------------------------------------------------------------------------------
% % -------------------------------------------------------------------------------
% %   --------------
% %   COARSE MESHES  !
% % ----------------
% [z,~,~,zNODES,Mbar,AREA_coarse_A] = CoarseMeshVariable(NameFileMeshCoarse,CentroidFA,COORrelA,nodesfA,COOR_FACE,AREA) ; 
% 
%  % Intersection of both sets 
% h = 1:length(DOFA); 
% h = setdiff(h,[z]) ; 
% 
% % Number of constraints 
% nconstr = 6  + length(z) ; 
% % number of unknowns 
% N = 2*length(z) ; 
% % Matrix of coefficients 
% A = sparse(nconstr,N) ; 
% b = zeros(nconstr,1) ; 
% 
% %%%%%% First row 
% % First column
% iini = 1;
% ifin = 6 ;
% jini = 1 ;
% jfin = length(z) ;
% A(iini:ifin,jini:jfin) =  R(z,:)'*Mbar ;
% b(iini:ifin) =  (R(z,:)'*Mbar)*(R(z,:)*a_A) ;
%  
%  
% 
% % IDENTITIES MATRICES 
% % --------------------
% iini = ifin+1; 
% ifin = nconstr ; 
% nrows = length(z) ; 
% jini = 1 ;
% jfin = nrows ;
% A(iini:ifin,jini:jfin) = speye(nrows,nrows) ;  
% 
% jini = jfin + 1 ; 
% jfin = jini + nrows -1 ; 
% A(iini:ifin,jini:jfin) = -speye(nrows,nrows) ;  
% b(iini:ifin) = R(z,:)*da ; 
% 
% DA = DOFA(z) ; 
% DB = DOFB(z) ; 
% D = [DA; DB] ; 
% 
% % %%% 
% % Abar = A(1:12,1:length(DOFA)) ; 
% % bA = b(1:12) ; 
% % bB = b(13:end) ; 
% 
% [~,r]=licols(A) ; %
% l = setdiff(1:2*length(z),r) ;
%  
% G = -A(:,r)\A(:,l) ; 
% uBAR = A(:,r)\b ; 
% DOFr = D(r) ; 
% DOFm = D(l) ; 
% 
% 
% 
%  
% 
%   
% 
% %  
% % 
% % R_z = R(z,:) ;
% % %Rast_z = Rast(z,:) ; 
% % % Rbar --->
% % Rbar_z = zeros(size(R_z)) ;
% % for idim =1:ndim
% %     INDLOC =idim:ndim:size(R_z,1) ;
% %     Rbar_z(INDLOC,:) = Mst_coarse*R_z(INDLOC,:) ;
% % end
% % 
% % 
% % % Select 6 linearly independent rows from Rbar_z
% % [~,r]=licols(Rbar_z') ; %
% % l = setdiff(1:length(z),r) ;
% % %
% % J = inv(Rbar_z(r,:)')*Rbar_z(l,:)' ;
% % b = inv(Rbar_z(r,:)')*Rbar_z'*R_z*a_A ;
% % DOFr = [DOFA(z(r)); DOFB(z(r)) ; DOFB(z(l)); DOFB(h)]  ;
% % DOFm = [DOFA(z(l)) ;DOFA(h) ];
% % G = sparse(length(DOFr),length(DOFm)) ;
% % uBAR = zeros(length(DOFr),1) ;
% % 
% % iini = 1;
% % ifin = size(J,1) ;
% % jini = 1 ;
% % jfin = size(J,2) ;
% % G(iini:ifin,jini:jfin) = - J ;
% % uBAR(iini:ifin) = b  ;
% % 
% % iini = ifin+1;
% % ifin = iini+size(J,1)-1 ;
% % G(iini:ifin,jini:jfin) = - J ;
% % uBAR(iini:ifin) = b - R_z(r,:)*da;
% % 
% % iini = ifin +1;
% % ifin = iini-1 + length(l) ;
% % G(iini:ifin,jini:jfin) = speye(length(l));
% % uBAR(iini:ifin) =   - R_z(l,:)*da;
% % 
% % iini = ifin +1;
% % jini = length(l) + 1;
% % G(iini:end,jini:end) = speye(length(h));
% % uBAR(iini:end) =   - R(h,:)*da;
% % 
