function [DOFr,dR] =  DirichletBCstandardRVE(NODES_LINES,ilines,ndim,DISPLOC) 


 % Old method (before 16-July-2019)
    DOFS = small2large(NODES_LINES{ilines},ndim) ; % Numbering of DOFs assuming ndim = ndimMAX
    nnodes = length(NODES_LINES{ilines}) ;
    % Known DOFs
    dR = [] ;
    DOFr = [] ;
    ndimLOC = length(DISPLOC) ;
    PRESCRIBED_ALL = 0 ;
    for idim = 1:ndim
        if idim <=ndimLOC
            if ~isempty(DISPLOC{idim})
                DOFr = [DOFr; (idim:ndim:length(DOFS))'] ;
                dR = [dR ; DISPLOC{idim}*ones(nnodes,1)] ;
                PRESCRIBED_ALL = 1; 
               % warning('Set it to 1')
            end
        elseif idim > ndimLOC
            if PRESCRIBED_ALL == 1
                DOFr = [DOFr; (idim:ndim:length(DOFS))'] ;
                dR = [dR ; 0*ones(nnodes,1)] ;
            end
        end
    end
    DOFr = DOFS(DOFr) ;