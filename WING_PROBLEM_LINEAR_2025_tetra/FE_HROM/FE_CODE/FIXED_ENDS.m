
rnod = cell(ndim,1) ; uPRES=cell(ndim,1)  ;

    TOL = ChooseTolerance(CN,COOR) ;

%
FIXED_END = INPUTS_LOC.FIXED_END ;
% Loop over FIXED_END
BoundaryNodes= unique(CONNECTb(:)) ;

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
    for idim = 1:ndim
        rnod{idim} =[ rnod{idim}  ; rnodBASE];
        uPRES{idim} =  [ uPRES{idim} ; zeros(size(rnodBASE))] ; %zeros(size(rnod{idim})) ;
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