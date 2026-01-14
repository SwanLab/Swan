% function[d, strainGLOgid, stressGLOgid,   React, posgp, CN, MaterialType, DATAOUT]  = ...
%     SolveElastFE_old(COOR,CN,TypeElement,TypeElementB, celasglo,  DOFr,dR,...
%     Tnod,CNb,fNOD,Fpnt,typePROBLEM,celasgloINV,DATA,CONNECTb,DOFm,Gb,MaterialType) ;
% 
% %%% This function returns the (nnode*ndim x 1) vector of nodal displacements (d),
% %%% as well as the arrays containing  the stresses (stressGLO) and strains (strainGLO) at all gauss
% %%% points
% % % INPUTS
% % --------------
% % 1. Finite element mesh
% % -------------------
% % COOR: Coordinate matrix (nnode x ndim)
% % CN: Connectivity matrix (nelem x nnodeE)
% % TypeElement: Type of finite element (quadrilateral,...)
% % TypeElementB: Type of boundary finite element (linear...)
% % -----------
% % 2. Material
% % -----------
% %  celasglo (nstrain x nstrain x nelem)  % Array of elasticity matrices
% %  celasgloINV (6 x 6 x nelem)  % Array of compliance matrices (3D)
% % -------------------------
% % 3. Dirichlet (Essential) Boundary Condition s
% % --------------------------------------------
% %  DOFr --> Set of Global Degrees of Freedom with prescribed displacements
% %  dR   --> Vector of prescribed displacements  (size(DOFr) = size(dR))
% % ---------------------------------------
% % 4. Neumann (natural) Boundary Conditions
% % -------------------------------------------
% % DISTRIBUTED LOADS
% % -------------------
% %  CNb: Cell array in which the cell CNb{idim} contains the connectivity matrix for the boundary elements
% %  of the traction boundaries in the idim direction
% %  Tnod: Cell array in which the entry Tnod{idim} features the vector with the prescribed traction at
% %   the nodes specified in CNb{idim}    (note that size(CNb{idim}) = size(Tnod{idim}))
% %  each cell of
% %  POINT LOADS
% % --------------------------------
% %  Fpnt  (nnode*ndime x 1):  Vector containing point forces applied on the
% %  nodes of the discretization
% % ----------------------------------
% % 5. Body force
% % ---------------
% %  fNOD: Vector containing the nodal values of the heat source function (nnode*ndime x1 )%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 6. Miscellaneous input DATA  -->
% % DATA.VECTcode = 1 --> VEctorized code
% %
% d=[]; strainGLO=[] ; stressGLO=[] ;posgp=[] ;
% % A) Global stiffness= matrix
% % ------------------------------
% disp('Computing stiffness matrix K ...')
% disp(['DATA.RECALCULATE_STIFFNESS=',num2str(DATA.RECALCULATE_STIFFNESS)]) ; 
% if isfield(DATA,'AREA')
%  AREA = DATA.AREA ; 
% else
%     AREA = [] ; 
% end
% if  DATA.RECALCULATE_STIFFNESS == 1
%     %   dbstop('53')
%     if DATA.VECTcode  == 0
%         % Standard way (elementwise)
%         K = ComputeK(COOR,CN,TypeElement, celasglo) ;
%         wSTs= [] ; wsT = [] ; XeALL = [] ;Cglo = [] ; Bst = [] ;
%         IndicesRenumberingElements = [] ; 
%     else
%         % Vectorized code
%         [K CN wSTs  XeALL Cglo Bst wST MaterialType IndicesRenumberingElements]=...
%             ComputeKvect(COOR,CN,TypeElement, celasglo,DATA,MaterialType) ;
%     end
%     %    dbstop('61')
%     if DATA.STORE_STIFFNESS ==1 |  DATA.STORE_STIFFNESS ==2
%         disp('Storing Stiffness Matrix...')
%         %     save(DATA.nameWORKSPACE,'K','CN','wSTs','XeALL','Cglo','Bst','wST','MaterialType','COOR','-append');
%         
%         if isfield(DATA,'densGLO')
%         density = DATA.densGLO ; 
%         else
%             density = [] ; 
%         end
%         save(DATA.nameWORKSPACE,'wST','Bst','wSTs','MaterialType','CN','Cglo','COOR','K','XeALL',...
%             'IndicesRenumberingElements','AREA','density','-append');
%         disp('Done')
%         
%         
%     end
%     
% else
%     
%     DATA = DefaultField(DATA,'nameWORKSPACE_Kstiff',DATA.nameWORKSPACE) ;
%     
%     disp('Retrieving Stiffness Matrix...')
%     %  dbstop('69')
%     load(DATA.nameWORKSPACE_Kstiff,'K','CN','wSTs','XeALL','Cglo','Bst','wST','MaterialType');
%     disp('Done')
% end
% 
% %save('K','K')
% 
% % -----------------------
% 
% %------------------------
% 
% % B) External force vector due to body forces
% % ------------------------------
% disp('Computing   external force vector due to body forces (Fb)...')
% if DATA.VECTcode  == 0
%     % Standard way (elementwise)
%     Fb = ComputeFb(COOR,CN,TypeElement, fNOD); Nst = [] ;
% else
%     % Vectorized code
%     % dbstop('90')
%     [Fb, Nst,  DATA, wSTs_RHS, posgp_RHS]= ComputeFbVect(COOR,CN,TypeElement, fNOD,DATA);
% end
% 
% if exist(DATA.nameWORKSPACE)==0
%     APPEND = '' ;
% else
%     APPEND = '-append' ;
% end
% 
% if DATA.STORE_STIFFNESS ==2 
%     save(DATA.nameWORKSPACE,'Nst','wSTs_RHS','fNOD','posgp_RHS','Fb',APPEND);
% elseif DATA.STORE_STIFFNESS ==1
%     save(DATA.nameWORKSPACE,'Fb',APPEND) ;
% end
% 
% 
% 
% % C)  External force vector due to   boundary tractions
% % ------------------------------
% disp('Computing  external force vector due to   boundary tractions ..')
% %dbstop('91')
% 
% 
% % error detected in the vectorized form
% % DATA_VECTcode =1 ;
% %dbstop('112')
% if DATA.VECTcode  == 0
%     % Standard way (elementwise)
%     Ftrac = FtracCOMP(COOR,CONNECTb,TypeElementB,Fpnt,Tnod);
% else
%     % Vectorized code
%     Ftrac = FtracCOMPvect(COOR,CNb,TypeElementB,Fpnt,Tnod,CONNECTb);
% end
% 
% if DATA.STORE_STIFFNESS ==1
%     save(DATA.nameWORKSPACE,'Ftrac','-append');
% elseif DATA.STORE_STIFFNESS ==2
%     save(DATA.nameWORKSPACE,'Ftrac','CNb','TypeElementB','Fpnt','Tnod','CONNECTb','-append')
% end
% 
% %dbstop('109')
% if DATA.CALCULATE_MASSMATRIX == 1
%     ndim = size(COOR,2) ;
%     
%     % if  DATA.RECALCULATE_STIFFNESS == 1
%     M= MassMatrix(DATA,Nst,wSTs_RHS,ndim) ;
%     % else
%     
%     disp('Storing Mass Matrix...')
%     save(DATA.nameWORKSPACE,'M','-append');
%     disp('Done')
%     %end
%     
% end
% 
% %save('M','M')
% 
% 
% % D) Solving for the vector of unknown displacements
% disp('Solving...')
% if DATA.NOCALCULATE_DISPLACEMENTS == 0
%     [d strainGLOgid stressGLOgid  React posgp DATAOUT] =...
%         SolveELAS(K,Fb,Ftrac,dR,DOFr,COOR,CN,TypeElement,celasglo,typePROBLEM,celasgloINV,...
%         Cglo,Bst,DOFm,Gb,DATA,wST,Nst) ;
% else
%     DATAOUT.stress = [] ;     strainGLOgid = [] ;     stressGLOgid = [] ;
%     d = [] ; React = [] ;
% end
% 
% % DATAOUT.strainGLO = strainGLO ;
% % DATAOUT.stressGLO = stressGLO ;
% DATAOUT.wSTs = wSTs ;
% DATAOUT.MaterialType = MaterialType ;
% DATAOUT.d = d ;
% %DATAOUT.Bst = Bst;
% stressGLO = DATAOUT.stress;
% DATA_INPUT_FE = DATA ; 
% save(DATA.nameWORKSPACE,'d','stressGLO','DOFr','TypeElement','posgp','DATA_INPUT_FE','-append')
% DATAOUT.nameWORKSPACE = DATA.nameWORKSPACE ;
