% See Implementation.pdf
% disp('BORRARRRRRRRRRRRRRR')
%load('tmp.mat')
%rnod = cell(ndim,1) ; uPRES=cell(ndim,1)  ;
% 

% 
TOL = ChooseTolerance(CN,COOR) ;
E = INPUTS_LOC.STRAIN;
NODESbound = unique(CONNECTb(:)); % All boundary nodes
COORbound = COOR(NODESbound,:) ;
xmin = min(COORbound(:,1)) ; xmin = xmin(1) ;
xmax = max(COORbound(:,1)) ; xmax = xmax(1) ;
ymin = min(COORbound(:,2)) ; ymin = ymin(1) ;
ymax = max(COORbound(:,2)) ; ymax = ymax(1) ;
zmax = [] ; zmin =[] ;
if size(COOR,2) ==2
    strain = [E(1) E(3); E(3) E(2) ] ; 
    
else
    strain = [E(1)  E(6)  E(5)
              E(6)  E(2)  E(4)
              E(5)  E(4)  E(3)];
          
          zmin = min(COORbound(:,3)) ; zmin = zmin(1) ;
zmax = max(COORbound(:,3)) ; zmax = zmax(1) ;
 
end

DATA = DefaultField(DATA,'TOL_deter_BOUNDARYNODES',[]) ; %.TOL_deter_BOUNDARYNODES
if isempty(DATA.TOL_deter_BOUNDARYNODES)
    TOL = ChooseTolerance(CN,COOR) ;
else
    TOL = DATA.TOL_deter_BOUNDARYNODES ;
end
 

[NODESln NormalPlanes]=  DetermineePlanesPeriodicHexaBCS(COORbound,TOL,xmax,xmin,ymax,ymin,zmax,zmin,NODESbound) ;
% ----------
BoundaryNodes = unique(cell2mat(NODESln)) ; 
COORbnd = COOR(BoundaryNodes,:) ; 

disp = strain*COORbnd' ; 
dR = disp(:) ; 
DOFr = small2large(BoundaryNodes,size(COOR,2)); 

% 
%  
% %dbstop('12')
% for iface =1:length(FIXED_END)
%     
%     FACE = FIXED_END{iface} ;
%     [dimREF signREF] = Faces2Index(FACE) ;
%     
%     
%     
%     % 1) List of nodes at which displacement is prescribed  in DIRECTION  i = 1
%     % All nodes pertaining to plane z = 0
%     % Boundary nodes
%     %  dbstop('42')
%     coorREF = min(signREF*COOR(BoundaryNodes,dimREF)) ; coorREF = signREF*coorREF(1) ;
%     rnodBASEloc = find(abs(COOR(BoundaryNodes,dimREF)-coorREF)<TOL) ;
%     rnodBASE = BoundaryNodes(rnodBASEloc) ;
%     
%     %%%% PRESCRIBED DISPLACEMENT 
%     %REFpoint = rnodBASE(1); % Reference point (ideally, this point should be the center of gravity)
%     
%     COOR_FACE = COOR(rnodBASE,:) ; 
%     COORref = sum(COOR_FACE,1)/size(COOR_FACE,1); % Center of gravity
%     
%     COORrel = bsxfun(@minus,COOR_FACE',COORref')'; % Relative coordinates
%     BasisUrb = ConstructBasisRigidBody(COORrel) ; % Basis matrix for rigid body motions
%     FACE = 1:size(BasisUrb,1) ;
%     uBARloc = PrescribedRIGIBbodyDisp(BasisUrb,FACE,DISP_RBODY{iface}) ; % Prescribed displacement
%     %%%%
%     
%     for idim = 1:ndim
%         rnod{idim} =[ rnod{idim}  ; rnodBASE];
%         uPRES{idim} =  [ uPRES{idim} ; uBARloc(idim:ndim:end)] ; %zeros(size(rnod{idim})) ;
%     end
% end
% 
% %dbstop('33')
% for i = 1:length(rnod)
%     [rnod{i}  iii jjj] = unique(rnod{i}  ) ;
%     uPRES{i} = uPRES{i}(iii) ;
% end
% 
% save(DATA.nameWORKSPACE,'rnod','-append')
% 
% 
% 
% DOFr = [] ; dR = [] ;
% for idim = 1:ndim
%     DOFr = [DOFr ; (rnod{idim}-1)*ndim+idim];
%     dR = [dR ; uPRES{idim}];
% end