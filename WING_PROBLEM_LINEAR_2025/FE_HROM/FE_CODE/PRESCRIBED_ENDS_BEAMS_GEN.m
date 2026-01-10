% clc
% clear all 
% warning('borra estoooooo')
% load('tmp.mat')
rnod = cell(ndim,1) ; uPRES=cell(ndim,1)  ;

a_A = INPUTS_LOC.DISPLACEMENTS_ENDS_RIGID_BODY{1}(:) ; % Rigid body amplitudes face 1
a_B = INPUTS_LOC.DISPLACEMENTS_ENDS_RIGID_BODY{2}(:) ; % Rigid body amplitudes face 2

if ndim == 2
    ROTind = [3] ;
else
    ROTind = [4:6] ; 
end

if any(cellfun(@isempty,a_A(ROTind)))
   error('Non-compatible boundary conditions ') 
end
if any(cellfun(@isempty,a_B(ROTind)))
   error('Non-compatible boundary conditions ') 
end

 
a = {a_A,a_B} ;  

%dbstop('12')
DATA = DefaultField(DATA,'angROTATION_FACE',[]) ;
DOMS = [1,size(CONNECTb,1)] ;
for iface =1:length(a)
    aLOC = a{iface} ; 
    INDZ = find(cellfun(@isempty,aLOC)) ;  
    aRIGB = aLOC ; 
    aRIGB(INDZ) = {0} ; 
    aRIGB = cell2mat(aRIGB) ; 
    % ROTATION
    % -------- 
    ROTMATRIX = LocalRotMatrix(iface,DATA,DOMS,ndim) ;  
    % 1) List of nodes at which displacement is prescribed  in DIRECTION  i = 1
    rnodBASE = unique(CONNECTb{DOMS(iface),iface}) ;
    
    %%%% PRESCRIBED DISPLACEMENT
    %REFpoint = rnodBASE(1); % Reference point (ideally, this point should be the center of gravity)
    
    COOR_FACE = COOR(rnodBASE,:) ;
    COORref = sum(COOR_FACE,1)/size(COOR_FACE,1); % Center of gravity    
    COORrel = bsxfun(@minus,COOR_FACE',COORref')'; % Relative coordinates    
    % ROTATION 
    % -------
    if ~isempty(ROTMATRIX)
        COORrel = (ROTMATRIX'*COORrel')' ; 
    end    
    BasisUrb = ConstructBasisRigidBody(COORrel) ; % Basis matrix for rigid body motions
    FACE = 1:size(BasisUrb,1) ;
    uBARloc = PrescribedRIGIBbodyDisp(BasisUrb,FACE,aRIGB) ; % Prescribed displacement
    %%%%    
    if iface==2 && ~isempty(DATA.angROTATION_FACE)        
        uBARloc = ROTMATRIX*reshape(uBARloc,ndim,[]) ;
        uBARloc = uBARloc(:) ;        
    end   
    
    for idim = 1:ndim
        if isempty(find(idim == INDZ))
        rnod{idim} =[ rnod{idim}  ; rnodBASE];
        uPRES{idim} =  [ uPRES{idim} ; uBARloc(idim:ndim:end)] ; %zeros(size(rnod{idim})) ;
        end
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