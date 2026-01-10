function DATA  = ReactionPlotONLINE_multiscale1D(DATA,OPERfe)
% Resultant reactions  % Change 2020-Jan-22
DATA = DefaultField(DATA,'REACTIONS_RESULTANTS_CALCULATE',[]) ;
DATA = DefaultField(DATA,'NODES_ENTITIES',[]) ;
DATA = DefaultField(DATA,'BasisUrb_ENTITIES',[]) ;

ReactionResultantDOFS = [] ;
DATA.ReactionResultants = [] ;

if ~isempty(DATA.REACTIONS_RESULTANTS_CALCULATE)   && ~isempty(DATA.NODES_ENTITIES)
    DATA.REACTIONS_RESULTANTS_CALCULATE = DefaultField(DATA.REACTIONS_RESULTANTS_CALCULATE,'ENTITY',1) ;
    DATA.REACTIONS_RESULTANTS_CALCULATE = ...
        DefaultField(DATA.REACTIONS_RESULTANTS_CALCULATE,'DIRECTION','x') ;
    
    ReactionResultantDOFS = cell(1,length(DATA.REACTIONS_RESULTANTS_CALCULATE));
    OPERATORreactions = cell(1,length(DATA.REACTIONS_RESULTANTS_CALCULATE));
    
    for ientityLOC = 1:length(DATA.REACTIONS_RESULTANTS_CALCULATE.ENTITY)
        ientity = DATA.REACTIONS_RESULTANTS_CALCULATE.ENTITY(ientityLOC) ;
        
        
        NODES = DATA.NODES_ENTITIES{ientity} ;
        DOFS = small2large(NODES,OPERfe.ndim) ;
        
        switch DATA.REACTIONS_RESULTANTS_CALCULATE.DIRECTION{ientityLOC}
            case {'x'}
                %   DOFSloc = DOFS(1:OPERfe.ndim:end) ;
                %  op
                OPERATORreactions{ientityLOC} = zeros(size(DOFS)) ; % 11-May-2020
                OPERATORreactions{ientityLOC}(1:OPERfe.ndim:end) = 1;
            case 'y'
                %   DOFSloc = DOFS(2:OPERfe.ndim:end) ;
                OPERATORreactions{ientityLOC} = zeros(size(DOFS)) ; % 11-May-2020
                OPERATORreactions{ientityLOC}(2:OPERfe.ndim:end) = 1;
                
            case {'Mz','M'}
             %   DOFSloc = DOFS  ;
                %    if ~isempty(DATA.BasisUrb_ENTITIES)
                if OPERfe.nstrain ==4 || OPERfe.nstrain ==3
                    % 2D problems
                    icomp = 3;
                else
                    icomp = 6;
                end
                OPERATORreactions{ientityLOC} = zeros(size(DOFS)) ; % 11-May-2020
                OPERATORreactions{ientityLOC}(icomp:OPERfe.ndim:end) = 1;
                %   else
                %      error('Rigid body modes not defined')
                %  end
                
            otherwise
                error('Option not implemented yet')
        end
        
        ReactionResultantDOFS{ientityLOC} = DOFS ;
        
    end
    
    DATA.ReactionResultantDOFS = ReactionResultantDOFS ;
    DATA.OPERATORreactions = OPERATORreactions;
    
    DATA.ReactionResultants = zeros(length(DATA.REACTIONS_RESULTANTS_CALCULATE.ENTITY),length(DATA.TIME_DISCRETIZATION));
    
    
end


