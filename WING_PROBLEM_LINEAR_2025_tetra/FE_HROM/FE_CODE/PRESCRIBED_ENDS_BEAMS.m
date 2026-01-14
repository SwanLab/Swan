
rnod = cell(ndim,1) ; uPRES=cell(ndim,1)  ;

a_A = INPUTS_LOC.DISPLACEMENTS_ENDS_RIGID_BODY{1}(:) ; % Rigid body amplitudes face 1
a_B = INPUTS_LOC.DISPLACEMENTS_ENDS_RIGID_BODY{2}(:) ; % Rigid body amplitudes face 2
if iscell(INPUTS_LOC.DISPLACEMENTS_ENDS_RIGID_BODY{1})
    INPUTS_LOC.DISPLACEMENTS_ENDS_RIGID_BODY{1} = cell2mat(INPUTS_LOC.DISPLACEMENTS_ENDS_RIGID_BODY{1}) ;
end
if iscell(INPUTS_LOC.DISPLACEMENTS_ENDS_RIGID_BODY{2})
    INPUTS_LOC.DISPLACEMENTS_ENDS_RIGID_BODY{2} = cell2mat(INPUTS_LOC.DISPLACEMENTS_ENDS_RIGID_BODY{2}) ;
end

%
% FUNinput.INPUTS.DISPLACEMENT_COND_FILE = 'PRESCRIBED_ENDS' ;
% FIXED_FACE = {'xMIN','xMAX'} ; %
% DISPLACEMENTS_ENDS_RIGID_BODY = {zeros(6,1),zeros(6,1)};
% DISPLACEMENTS_ENDS_RIGID_BODY{2}(6) =0.1; % Rotation around z-axis
% FUNinput.INPUTS.FIXED_END =FIXED_FACE ;
% FUNinput.INPUTS.DISPLACEMENTS_ENDS_RIGID_BODY =DISPLACEMENTS_ENDS_RIGID_BODY ;
%
%
%TOL = ChooseTolerance(CN,COOR) ;
INPUTS_LOC = DefaultField(INPUTS_LOC,'FIXED_END',[1 2]) ; 
FIXED_END = INPUTS_LOC.FIXED_END ; 
if isempty( INPUTS_LOC.DISPLACEMENTS_ENDS_RIGID_BODY{1}) 
    FIXED_END = 2 ; 
elseif isempty( INPUTS_LOC.DISPLACEMENTS_ENDS_RIGID_BODY{2})
     FIXED_END = 1 ;
elseif isempty( INPUTS_LOC.DISPLACEMENTS_ENDS_RIGID_BODY{2}) && isempty( INPUTS_LOC.DISPLACEMENTS_ENDS_RIGID_BODY{1})
    error('Statically undetermined problem')
end


DISP_RBODY = INPUTS_LOC.DISPLACEMENTS_ENDS_RIGID_BODY ;
% Loop over FIXED_END
%BoundaryNodes= unique(CONNECTb(:)) ;

%BasisRrb = ConstructBasisRigidBody(COORrel) ;
%        uBARloc = PrescribedRIGIBbodyDisp(BasisUrb,FACE,uBAR_0)  ;

%dbstop('12')
DATA = DefaultField(DATA,'angROTATION_FACE',[]) ;
DOMS = [1,size(CONNECTb,1)] ;

DOFA = [] ; 
DOFB = [] ; 

RotMatrixA = [] ; 
RotMatrixB = [] ; 

NODES_ENTITIES = cell(size(FIXED_END)) ; 
BasisUrb_ENTITIES = cell(size(FIXED_END)) ; 
for iface =1:length(FIXED_END)
    
        % ROTATION
    % -------- 
    ROTMATRIX = LocalRotMatrix(iface,DATA,DOMS,ndim) ;  
    % 1) List of nodes at which displacement is prescribed  in DIRECTION  i = 1
    rnodBASE = unique(CONNECTb{DOMS(iface),FIXED_END(iface)}) ;
    NODES_ENTITIES{iface} = rnodBASE ; 
 
    
    %%%% PRESCRIBED DISPLACEMENT
    %REFpoint = rnodBASE(1); % Reference point (ideally, this point should be the center of gravity)
    
    
    
    COOR_FACE = COOR(rnodBASE,:) ;
    
     %  
 
     compute_exact_centroid =1 ; 
     
     if compute_exact_centroid == 0
         % Before 10-June-2022
    COORref = sum(COOR_FACE,1)/size(COOR_FACE,1); % Center of gravity    
     else
       %  error('REVISA ESTO !!! ')
       
          [COORref,AREA,Mst] =CentroidGeometricMassMatrixNEW(COOR,rnodBASE,CONNECTb{DOMS(iface),FIXED_END(iface)},TypeElementB) ; 
     end
    
    
    
    COORrel = bsxfun(@minus,COOR_FACE',COORref')'; % Relative coordinates    
    % ROTATION 
    % -------
    if ~isempty(ROTMATRIX)
        COORrel = (ROTMATRIX'*COORrel')' ; 
    end    
    BasisUrb = ConstructBasisRigidBody(COORrel) ; % Basis matrix for rigid body motions
    BasisUrb_ENTITIES{iface} = BasisUrb ; 
    
       if iface ==1 
        DOFA = small2large(rnodBASE,ndim) ; 
        BasisUrbA = BasisUrb ; 
        RotMatrixA = ROTMATRIX ; 
    else
        DOFB = small2large(rnodBASE,ndim) ; 
          BasisUrbB = BasisUrb ; 
             RotMatrixB = ROTMATRIX ; 
    end
    
    FACE = 1:size(BasisUrb,1) ;
    uBARloc = PrescribedRIGIBbodyDisp(BasisUrb,FACE,DISP_RBODY{iface}) ; % Prescribed displacement
    %%%%    
    if iface==2 && ~isempty(DATA.angROTATION_FACE)        
        uBARloc = ROTMATRIX*reshape(uBARloc,ndim,[]) ;
        uBARloc = uBARloc(:) ;        
    end   
    
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

%save(DATA.nameWORKSPACE,'rnod','-append')



DOFr = [] ; dR = [] ;
for idim = 1:ndim
    DOFr = [DOFr ; (rnod{idim}-1)*ndim+idim];
    dR = [dR ; uPRES{idim}];
end