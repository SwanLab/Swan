% dbstop('2')
% load('tmp.mat')
% disp('BORRARRRRRRRRRRRRRR')
rnod = cell(ndim,1) ; uPRES=cell(ndim,1)  ;
% 
% FUNinput.INPUTS.DISPLACEMENT_COND_FILE = 'PRESCRIBED_ENDS' ; 
% FIXED_FACE = {'xMIN','xMAX'} ; %
% DISPLACEMENTS_ENDS_RIGID_BODY = {zeros(6,1),zeros(6,1)}; 
% DISPLACEMENTS_ENDS_RIGID_BODY{2}(6) =0.1; % Rotation around z-axis
% FUNinput.INPUTS.FIXED_END =FIXED_FACE ; 
% FUNinput.INPUTS.DISPLACEMENTS_ENDS_RIGID_BODY =DISPLACEMENTS_ENDS_RIGID_BODY ; 
% 
% 
TOL = ChooseTolerance(CN,COOR) ;
FIXED_END = INPUTS_LOC.FIXED_END ;
DISP_RBODY = INPUTS_LOC.DISPLACEMENTS_ENDS_RIGID_BODY ; 
% Loop over FIXED_END
BoundaryNodes= unique(CONNECTb(:)) ;

if isempty(BoundaryNodes)
   error('Boundary mesh information missing... You should generate the mesh boundary before exporting the .msh file') 
end

%BasisRrb = ConstructBasisRigidBody(COORrel) ; 
%        uBARloc = PrescribedRIGIBbodyDisp(BasisUrb,FACE,uBAR_0)  ; 

%dbstop('12')
for iface =1:length(FIXED_END)
    
    FACE = FIXED_END{iface} ;
    [dimREF signREF] = Faces2Index(FACE) ;
    
    
    
    % 1) List of nodes at which displacement is prescribed  in DIRECTION  i = 1
    % All nodes pertaining to plane z = 0
    % Boundary nodes
    %  dbstop('42')
    coorREF = min(signREF*COOR(BoundaryNodes,dimREF)) ; coorREF = signREF*coorREF(1) ;
    rnodBASEloc = find(abs(COOR(BoundaryNodes,dimREF)-coorREF)<TOL) ;
    rnodBASE = BoundaryNodes(rnodBASEloc) ;
    
    %%%% PRESCRIBED DISPLACEMENT 
    %REFpoint = rnodBASE(1); % Reference point (ideally, this point should be the center of gravity)
    
    COOR_FACE = COOR(rnodBASE,:) ; 
    COORref = sum(COOR_FACE,1)/size(COOR_FACE,1); % Center of gravity
    
    COORrel = bsxfun(@minus,COOR_FACE',COORref')'; % Relative coordinates
    BasisUrb = ConstructBasisRigidBody(COORrel) ; % Basis matrix for rigid body motions
    FACE = 1:size(BasisUrb,1) ;
    uBARloc = PrescribedRIGIBbodyDisp(BasisUrb,FACE,DISP_RBODY{iface}) ; % Prescribed displacement
    %%%%
    
    for idim = 1:ndim
        rnod{idim} =[ rnod{idim}  ; rnodBASE];
        uPRES{idim} =  [ uPRES{idim} ; uBARloc(idim:ndim:end)] ; %zeros(size(rnod{idim})) ;
    end
end

%dbstop('33')
for i = 1:length(rnod)
    [rnod{i}  iii jjj] = unique(rnod{i}  ) ;
    uPRES{i} = uPRES{i}(iii) ;
end

save(DATA.nameWORKSPACE,'rnod','-append')



DOFr = [] ; dR = [] ;
for idim = 1:ndim
    DOFr = [DOFr ; (rnod{idim}-1)*ndim+idim];
    dR = [dR ; uPRES{idim}];
end