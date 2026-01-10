% function [COOR,CN,TypeElement,CONNECTb,TypeElementB,Materials,MATERIALNEW,DOMAINVAR]= ...
%     MeshGenerationRepeat3D_2repeat(NAMEMSH,DATA,MATERIAL)
% % Automatic generation of meshes by tiling copies (translation in two
% % direction). See DomainDecom_SVD.m
% % Copy of MeshGenerationRepeat3D MeshGenerationRepeat3D_faces ,  (July  2018)
% % JaHO
% %
% % OUTPUTS:
% % --------
% %
% % DOMAINVAR.ListElements = ListElementsDom ; % List of elements of each domain
% % DOMAINVAR.NODES_faces12 = NODES_faces12 ; % List of connecting nodes of each domain (face 1  and face 2), so that
% % the (y,z) coordinates of NODES_faces12{idom,iface} and
% % NODES_faces12{jdom,iface} are the same for any "idom" and any "iface"
% % ----------------------------------------------------------------
% %addpath('DATA_input')
% % INPUTS
% % ---------
% if nargin ==0
%     load('tmp3.mat')
% end
% %%%%
% 
% nDOM = DATA.MakeMeshByRepetition.nDOM ; % number of domains per direction
% 
% %% nDOM(1) copies along the x-direction, nDOM(2) along the y-direction
% % ------------------------------------------------------------------------------
% disp('REpeating cells ...')
% disp('----------------------------')
% RVE.NAME =NAMEMSH ;
% DATA3D = GeometryRVE(RVE,DATA) ; % Reading mesh and node faces (already paired for faces 1, 2. and 3,4  )
% 
% DATA = DefaultField(DATA,'NAME_WS_MESH_GIVEN',[]) ;
% 
% %
% COOR = DATA3D.COOR ; % Coordinates
% CN = DATA3D.CN ; % Connectivities
% TypeElement = DATA3D.TypeElement   ;
% TypeElementB = DATA3D.TypeElementB ; % Type boundary element
% if ~isempty(DATA3D.MaterialType)
%     MaterialType = DATA3D.MaterialType ;
% else
%     MaterialType = ones(size(CN,1),1) ;
% end
% CONNECTb = DATA3D.CNb ;  % Face connectivities
% 
% MATERIALNEW.PLY = MATERIAL.PLY ;  % Set of materials defined in the input data
% 
% nmat = length(unique(MaterialType)) ;
% if nmat ~=length(MATERIAL.PLY)
%     error('The number of materials should coincide with the number of materials defined in the input file')
% end
% NODESfaces = DATA3D.NODES_FACES ;  % All face nodes (labels from GID)
% 
% % We only take into account  boundary connectivities of the
% % faces specified by NODESfaces (by the user, in GID)
% CONNECTb_faces = cell(1,length(NODESfaces)) ;
% for iface = 1:length(NODESfaces)
%     [dummy, setBelemLOC]= ElemBnd(CONNECTb,NODESfaces{iface}); % elements face "iface"
%     % Connectivities faces f1 and f2
%     CONNECTb_faces{iface} = CONNECTb(setBelemLOC,:) ;
%     % The above is the connectivity matrix for the nodes of face "iface"
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%%% Repetitition along the x-direction
% % -------------------------------------
% nDOMglo = nDOM ;
% 
% nDOM = nDOMglo(1) ;
% iface1 = 1;
% iface2  = 3 ;
% DIRECTION = 1;
% NODESf1 = DATA3D.NODES_FACES{iface1} ;
% NODESf2 = DATA3D.NODES_FACES{iface2} ;
% [COOR_x,CN_x,CONNECTb_x,MaterialType_x,MATERIALNEW_x,DOMAINVAR_x] = RepMeshALong1Dimension(nDOM,MATERIAL,MaterialType,CONNECTb_faces,CN,COOR,...
%     DATA3D,iface1,iface2,DIRECTION,nmat,NODESf1,NODESf2,NODESfaces) ;
% IMPRIMIR = 0; 
% if IMPRIMIR == 1
%  NameFile_msh = '/home/joaquin/Desktop/CURRENT_TASKS/COMPOSITE_MATERIALS_DOCENCIA/APUNTES_LATEX/DOCUMENTOS_anexos/MATLAB/ELASTOSTATIC_GEN/prueba.msh'
%  IND_ELEM_MESHES = GidMesh2DFE_multi(NameFile_msh,COOR_x,{CN_x},'',{MaterialType_x},{TypeElement},{'Prueba'});
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Recovering actual order of nodes face 2 and 4.
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NODESf2 = 2 ; 
% NODESf4  =4 ; 
% NODES_faces_24 = RecoveringListNodes24(DATA3D,CONNECTb_x,NODESf2,NODESf4) ;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Now we have to repeat the new mesh along the y-direction
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % We are going to use again "RepMeshALong1Dimension", but changing the
% % inputs
% nDOM = nDOMglo(2) ;
% iface1 = 2;
% iface2  = 4 ;
% DIRECTION = 2;
% 
% NODESf2 = cell2mat(NODES_faces_24(:,1));
% NODESf4 = cell2mat(NODES_faces_24(:,2)) ;
% 
% CONNECTb = cell(1,size(CONNECTb_x,2)) ;
% for iface = 1:size(CONNECTb,2)
%     CONNECTb{iface} = cell2mat(CONNECTb_x(:,iface)) ;
% end
% 
% 
% 
% if nDOM >1
%     nmat = length(unique(MaterialType_x)) ;
%     
%     [COOR,CN,CONNECTb_y,Materials,MATERIALNEW,DOMAINVAR_y] =...
%         RepMeshALong1Dimension(nDOM,MATERIALNEW_x,MaterialType_x,CONNECTb,CN_x,COOR_x,...
%         DATA3D,iface1,iface2,DIRECTION,nmat,NODESf2,NODESf4) ;
%     
% %     NameFile_msh = '/home/joaquin/Desktop/CURRENT_TASKS/COMPOSITE_MATERIALS_DOCENCIA/APUNTES_LATEX/DOCUMENTOS_anexos/MATLAB/ELASTOSTATIC_GEN/prueba.msh'
% %     IND_ELEM_MESHES = GidMesh2DFE_multi(NameFile_msh,COOR,{CN},'',{Materials},{TypeElement},{'Prueba'});
%     
%     % Now we have to build up the CONNECT_b and DOMAINVAR arrays.
%     % Let us start by the CONNECTb array.  (ndomx x ndomy x nfaces )
%     CONNECTb = RecoverCONNECTb(CONNECTb_y,nDOMglo) ;
%     
%     %% Next we do the same with the nodes of faces 1,2, 3 and 4
%     % Let us start by the nodes of faces 1 and 3
%     % We have the list corresponding to the first x-row --->
%     indNODESf1= 1 ; 
%     indNODESf3  =3 ; 
%     NODES_face13_1 =  DOMAINVAR_x.NODES_faces12;
%     NODES_faces_13 = RecoveringListNodes13(NODES_face13_1,CONNECTb,indNODESf1,indNODESf3) ;
%     % Now for faces 3 and 4. In his case we have at our disposal
%     % DOMAINVAR_y
%     NODES_faces24 = DOMAINVAR_y.NODES_faces12 ;
%     NODES_faces24 = RecoveringListNodes24_2(NODES_faces24,nDOMglo) ;
%     % New variable: NODES_faces  (for  faces 1,2,3,4)
%     NODES_faces = cell(nDOMglo(1),nDOMglo(2),4) ;
%     NODES_faces(:,:,[1,3]) = NODES_faces_13 ;
%     NODES_faces(:,:,[2,4]) = NODES_faces24 ;
%     DOMAINVAR.NODES_faces = NODES_faces ;
%     
%     %%% Now list of elements and nodes
%     % -------------------------------------------
%     % We have:
%     %DOMAINVAR_x.ListElements
%     %DOMAINVAR_y.ListElements
%     DOMAINVAR.ListElements = RecoveringListElements(DOMAINVAR_x.ListElements ,DOMAINVAR_y.ListElements)  ;
%     DOMAINVAR.ListNodesDom = RecoveringListElements(DOMAINVAR_x.ListNodesDom ,DOMAINVAR_y.ListNodesDom)  ;
%     
%     
% else
%     COOR = COOR_x ;
%     CN = CN_x ;
%      CONNECTb = cell(nDOMglo(1),nDOMglo(2),size(CONNECTb_x,2)) ; 
%      
%      for iface= 1:size(CONNECTb_x,2)
%      CONNECTb(:,1,iface) = CONNECTb_x(:,iface) ; 
%      end
%         
%      
%     Materials = MaterialType_x;
%     MATERIALNEW = MATERIALNEW_x
%     DOMAINVAR = DOMAINVAR_x ;
%     
%     NODES_faces = cell(nDOMglo(1),nDOMglo(2),4) ;
%     NODES_faces(:,:,[1 3]) = DOMAINVAR_x.NODES_faces13 ;
%     NODES_faces(:,:,[2 4]) = NODES_faces_24 ;
%     DOMAINVAR.NODES_faces = NODES_faces ;
%     
%     
% end
% 
% IMPRIMIR = 1 ; 
% if IMPRIMIR == 1
%  NameFile_msh = '/home/joaquin/Desktop/CURRENT_TASKS/COMPOSITE_MATERIALS_DOCENCIA/APUNTES_LATEX/DOCUMENTOS_anexos/MATLAB/ELASTOSTATIC_GEN/prueba.msh'
%  IND_ELEM_MESHES = GidMesh2DFE_multi(NameFile_msh,COOR,{CN},'',{Materials},{TypeElement},{'Prueba'});
% end
% 
% end
% 
% 
% function [COOR,CN,CONNECTb,Materials,MATERIALNEW,DOMAINVAR] = RepMeshALong1Dimension(nDOM,MATERIAL,MaterialType,CONNECTb_faces,CN,COOR,...
%     DATA3D,iface1,iface2,DIRECTION,nmat,NODES_f1,NODES_f2,NODESfaces)
% 
% NODES_faces12 = cell(nDOM,2) ;
% 
% 
% 
% 
% if nDOM== 1
%     % There is just one single domain
%     % --------------------------------
%     MATERIALNEW =    MATERIAL   ;
%     Materials =  MaterialType   ;
%     CONNECTb = cell(nDOM,length(CONNECTb_faces))  ;
%     for iface = 1:length(NODESfaces)
%         CONNECTb{1,iface} = CONNECTb_faces{iface} ; % Each entry should be a cell itself
%     end
%     ListElementsDom = {1:size(CN,1)} ;
%     ListNodesDom = {1:size(COOR,1)} ;
%     NODES_faces12{1,1} = NODES_f1 ;
%     NODES_faces12{1,2} = NODES_f2 ;
%     
%     
% %     idir = 1;
% %     COOR(:,idir) = DATA3D.CENTRf{iface1}(idir) + ( COOR(:,idir) - DATA3D.CENTRf{iface1}(idir)) ;
% %     
%     
% else
%     % -----------
%     f1NOD = NODES_f1;  % Nodes face 1, domain 1
%     f2NOD = NODES_f2;  % Nodes face 2, domain 1,
%     % Notice that f1NOD(i) matches with f2NOD(i)
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Global arrays
%     COORglo = cell(nDOM,1) ;  % Global coordinates
%     MATERIALNEW.PLY = MATERIAL.PLY ; 
%     
%     % Translation vector
%     translationVECTOR = COOR(f2NOD(1),:)-COOR(f1NOD(1),:);  % This should be changed for curved elements
%     
%     COOR1stdomain = COOR ;
%     idir = DIRECTION;
%     %%%% Initial point for the next domain
%     xINI = DATA3D.CENTRf{iface1} ;
%     xINI(idir) =  DATA3D.CENTRf{iface1}(idir) +   translationVECTOR(idir) ;
%     
%     COORglo{1} = COOR1stdomain ;
%     CNglo   = cell(nDOM,1) ;  % Global connectivities, all nodes
%     CNglo{1} = [CN]  ;
%     CONNECTbGLO=  cell(nDOM,length(CONNECTb_faces)) ; % Global connectivities, boundary nodes
%     CONNECTbGLO(1,:) =  CONNECTb_faces;
%     
%     % We create also an array for storing the nodes of faces F1 and F2 of
%     % each domain.
%     NODES_faces12 = cell(nDOM,2) ;
%     NODES_faces12{1,1} = f1NOD ;
%     NODES_faces12{1,2} = f2NOD ;
%     
%     % ---------------------------
%     Materials =cell(nDOM,1) ;
%     Materials{1} =MaterialType;
%     % --------------------------
%     ListElementsDom = cell(nDOM,1) ;
%     ListElementsDom{1} = 1:size(CN,1) ;  % List of elements domain 1
%     ListNodesDom = cell(nDOM,1) ;
%     ListNodesDom{1} = 1:size(COOR,1) ;  % List of nodes domain 1
%     
%     %% Loop over domains
%     nnodes = size(COOR,1) ;
%     ndim = size(COOR,2) ;
%     
%     %%% Coordinate relative to centroid face 1
%     COOR_REL = COOR ; %
%     for idim = 1:size(COOR,2)
%         COOR_REL(:,idim) =    COOR(:,idim) - DATA3D.CENTRf{iface1}(idim) ;
%     end
%     
%     for e = 2:nDOM
%         disp(['Cell =',num2str(e)])
%         
%         
%         %  -----------------------
%         % Coordinates domain e ***
%         % ------------------------
%         % Computation of the global coordinates of this domain
%         COORglo{e} =  zeros(size(COOR))  ;
%         idir = DIRECTION;
%         % COORglo{e}(:,idir)=  xINI(idir) + COOR_REL(:,idir) ;
%         for idim = 1:ndim
%             COORglo{e}(:,idim)=  xINI(idim) + COOR_REL(:,idim) ;
%         end
%         % Updating the reference point
%         xINI(idir) =  xINI(idir)  +   (translationVECTOR(idir)) ;
%         
%         
%         
%         %% List of nodes
%         % -----------------------------------------------------
%         ListNodesDom{e} = ListNodesDom{e-1}+nnodes ; % List of nodes domain ( for ensuring same numbering
%         % of nodes for all domains)
%         % Connectivity matrix domain. All elements
%         %----------------------------------------------
%         CNnew = CNglo{e-1} + nnodes; % We sum up the number of nodes of each domain
%         ListElementsDom{e} = ListElementsDom{e-1}+size(CNnew,1) ;
%         % Boundary elements (for each labeled face)
%         CNbLOC = CONNECTbGLO(e-1,:) ;
%         for iface =1:length(CONNECTb_faces)
%             CNbLOC{iface} = CONNECTbGLO{e-1,iface} + nnodes;
%         end
%         % -------------------------------------------
%         %         f1NEW = f1NOD + (e-1)*nnodes ;  % Numbering of new face 1
%         %         f2OLD = f2NOD + (e-2)*nnodes ;  % Numbering of old face 2
%         f1NEW = NODES_faces12{e-1,1} + nnodes;
%         f2OLD = NODES_faces12{e-1,2} ;
%         % -------------------------------------------------------------
%         NODES_faces12{e,1} = f2OLD ;   % Face 1,
%         NODES_faces12{e,2} = NODES_faces12{e-1,2} + nnodes ;   % Face 3
%         % ----------
%         
%         % Renumbering: we have to replace nodes (f1NEW) by f2NOD
%         % in CNnew and CNbLOC
%         ListNodesDom{e}(f1NOD) =  ListNodesDom{e-1}(f2NOD) ;
%         
%         CNnewREN  = CNnew ;
%         CNbLOCnew = CNbLOC ;
%         for ifacen = 1:length(f1NOD)   % This may be improved (make it more efficient)...
%             nodeLOC = f1NEW(ifacen);
%             % Replacing it in the connectivity matrix
%             INDnodes = find(CNnew==nodeLOC) ;
%             CNnewREN(INDnodes) = f2OLD(ifacen) ;
%             % Replacing it in the face connectivity matrix
%             for iface = 1:length(CNbLOCnew)
%                 INDnodesb = find(CNbLOC{iface}==nodeLOC) ;
%                 CNbLOCnew{iface}(INDnodesb) =f2OLD(ifacen) ;
%             end
%         end
%         
%         CNglo{e} =  CNnewREN ;
%         CONNECTbGLO(e,:) = CNbLOCnew ; %
%         Materials{e} = Materials{e-1} + nmat ;
%         for imat = 1:length(MATERIAL.PLY)
%             matREF = (e-1)*nmat ;
%             MATERIALNEW.PLY(matREF+imat) =  MATERIAL.PLY(imat) ;
%         end
%         
%     end
%     
%     Materials = cell2mat(Materials) ;
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     disp('REnumbering...')
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Numeration of nodes is not consecutive. We turn it consecutive in
%     % what follows
%     CNglo = cell2mat(CNglo) ;
%     NODES =  unique(CNglo);  % list of "old" node numbers, sorted from small to large
%     COORglo= cell2mat(COORglo) ;
%     COOR = COORglo(NODES,:) ;          % Global matrix of coordinates
%     NODES_new = 1:length(NODES) ;      % List of "new" node numbers, from 1 to length(NODES)
%     % Interior CN
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     CN = RenumberConnectivities(CNglo,NODES_new) ;
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Boundary CN
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     CONNECTb = CONNECTbGLO ;
%     for idom = 1:size(CONNECTbGLO,1)
%         for iface = 1:size(CONNECTbGLO,2)
%             NODESbnd = unique(CONNECTbGLO{idom,iface}) ;
%             [~,NODES_bnd,~] = intersect(NODES,NODESbnd) ;
%             CONNECTb{idom,iface}= RenumberConnectivities(CONNECTbGLO{idom,iface},NODES_bnd) ;
%         end
%         
%     end
%     
%     %% List of nodes FACES 1 and 2
%     for idom = 1:nDOM
%         for iface =1:size(NODES_faces12,2)
%             [~,NODES_bnd,~] = intersect(NODES,NODES_faces12{idom,iface}) ;
%             NODES_faces12{idom,iface}= RenumberConnectivities(NODES_faces12{idom,iface},NODES_bnd) ;
%         end
%     end
%     
%     %% List of all nodes
%     for idom = 1:nDOM
%         [~,NODES_dom,~] = intersect(NODES,ListNodesDom{idom}) ;
%         ListNodesDom{idom}= RenumberConnectivities(ListNodesDom{idom},NODES_dom) ;
%     end
%     
% end
% 
% DOMAINVAR.ListElements = ListElementsDom ; % List of elements of each domain
% DOMAINVAR.NODES_faces12 = NODES_faces12 ; % List of connecting nodes of each domain (face 1  and face 2), so that
% % the (y,z) coordinates of NODES_faces12{idom,iface} and
% % NODES_faces12{jdom,iface} are the same for any "idom" and any "iface"
% DOMAINVAR.ListNodesDom = ListNodesDom ;
% 
% 
% end
% 
% 
% 
% 
% 
% 
% 
% function  NODES_faces_24 = RecoveringListNodes24(DATA3D,CONNECTb,NODESf2,NODESf4)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Recovering actual order of nodes face 3 and 4.
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f2 = DATA3D.NODES_FACES{NODESf2} ; % Original mesh --> Reference face (face 3).
% f4 = DATA3D.NODES_FACES{NODESf4} ;
% % The above lists of nodes are sorted so that
% % nodes f3(i) and f4(i) are paired.
% % However, in the above function RepMeshALong1Dimension, the output is the
% % set of boundary connectivities. We have the lists of all nodes of faces 3
% % and 4 for all domains. However, we do not have at our disposal the actual
% % order in which these nodes have to be placed to be "paired"
% % -------------------
% idom = 1 ;
% IND_NODES_f2 = zeros(size(f2)) ;
% IND_NODES_f4= zeros(size(f4)) ;
% 
% for inode = 1:length(f2)
%     ORDERNODE =  find(CONNECTb{idom,NODESf2} == f2(inode)) ;
%     IND_NODES_f2(inode) = ORDERNODE(1) ;
%     ORDERNODE =  find(CONNECTb{idom,NODESf4} == f4(inode)) ;
%     IND_NODES_f4(inode) = ORDERNODE(1) ;
% end
% 
% % Therefore
% nDOM = size(CONNECTb,1) ;
% NODES_faces_24 = cell(nDOM,2) ;
% for idom = 1:nDOM
%     LISTNODES = CONNECTb{idom,NODESf2}(:) ;
%     NODES_faces_24{idom,1} = LISTNODES(IND_NODES_f4) ;
%     LISTNODES = CONNECTb{idom,NODESf4}(:) ;
%     
%     NODES_faces_24{idom,2} = LISTNODES(IND_NODES_f4) ;
% end
% end
% 
% function  CONNECTb = RecoverCONNECTb(CONNECTb_y,nDOMglo)
% 
% nfaces =size(CONNECTb_y,2) ;
% CONNECTb = cell(nDOMglo(1),nDOMglo(2),nfaces) ;
% % Establishing equivalence ---
% for idomy = 1:nDOMglo(2) % Loop over domains direction y
%     CNdomy = CONNECTb_y(idomy,:) ;
%     for iface = 1:nfaces
%         CNface = CNdomy{iface} ;
%         nrows = size(CNface,1)/nDOMglo(1) ;
%         nrows = nrows*ones(nDOMglo(1),1) ;
%         ncols = size(CNface,2) ;
%         NEWCN = mat2cell(CNface,nrows,ncols) ;
%         CONNECTb(:,idomy,iface) = NEWCN ;
%     end
%     
% end
% 
% end
% 
% function      NODES_faces_12 = RecoveringListNodes13(NODES_face12_1,CONNECTb,indNODESf1,indNODESf2)
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f1 = NODES_face12_1{1,1} ; % Original mesh --> Reference face (face 3).
% f2= NODES_face12_1{1,2};
% 
% idomx = 1 ;
% idomy = 1;
% IND_NODES_f1 = zeros(size(f1)) ;
% IND_NODES_f2= zeros(size(f2)) ;
% 
% for inode = 1:length(f1)
%     ORDERNODE =  find(CONNECTb{idomx,idomy,indNODESf1} == f1(inode)) ;
%     IND_NODES_f1(inode) = ORDERNODE(1) ;
%     ORDERNODE =  find(CONNECTb{idomx,idomy,indNODESf2} == f2(inode)) ;
%     IND_NODES_f2(inode) = ORDERNODE(1) ;
% end
% 
% % Therefore
% nDOMx = size(CONNECTb,1) ;
% nDOMy = size(CONNECTb,2) ;
% NODES_faces_12 = cell(nDOMx,nDOMy,2) ;
% for idom = 1:nDOMx
%     for jdom = 1:nDOMy
%         LISTNODES = CONNECTb{idom,jdom,indNODESf1}(:) ;
%         NODES_faces_12{idom,jdom,1} = LISTNODES(IND_NODES_f1) ;
%         LISTNODES = CONNECTb{idom,jdom,indNODESf2}(:) ;
%         
%         NODES_faces_12{idom,jdom,2} = LISTNODES(IND_NODES_f2) ;
%     end
% end
% 
% 
% end
% 
% function  NODES_faces24 = RecoveringListNodes24_2(NODES_faces24_inp,nDOMglo) ;
% 
% nfaces =size(NODES_faces24_inp,2) ;
% NODES_faces24 = cell(nDOMglo(1),nDOMglo(2),nfaces) ;
% % Establishing equivalence ---
% for idomy = 1:nDOMglo(2) % Loop over domains direction y
%     NODES = NODES_faces24_inp(idomy,:) ;
%     for iface = 1:nfaces
%         NODESloc = NODES{iface} ;
%         nrows = size(NODESloc,1)/nDOMglo(1) ;
%         nrows = nrows*ones(nDOMglo(1),1) ;
%         ncols = size(NODESloc,2) ;
%         NEWCN = mat2cell(NODESloc,nrows,ncols) ;
%         NODES_faces24(:,idomy,iface) = NEWCN' ;
%     end
%     
% end
% 
% 
% end
% 
% 
% function      ListElements = RecoveringListElements(ListElementsX,ListElementsY)
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%
% ndomx = length(ListElementsX) ;
% ndomy = length(ListElementsY) ;
% ListElements = cell(ndomx,ndomy) ;
% INDS_DOMX = cell(ndomx,1) ;
% idomy = 1;
% for idomx = 1:ndomx
%     ELEMS = ListElementsX{idomx} ;
%     INDS_DOMX{idomx} = zeros(size(ELEMS)) ;
%     
%     for ielem = 1:length(ELEMS)
%         ORDERELEM =  find(ListElementsY{idomy} == ELEMS(ielem)) ;
%         INDS_DOMX{idomx} (ielem) = ORDERELEM(1) ;
%     end
% end
% 
% % Therefore
% for idomx = 1:ndomx
%     for idomy = 1:ndomy
%         ListElements{idomx,idomy} = ListElementsY{idomy}(INDS_DOMX{idomx}) ;
%     end
% end
% 
% 
% 
% end
