function DATA  = ReactionPlotONLINE_multiscale(DATA,OPERfe)
% Resultant reactions  % Change 2020-Jan-2t
if nargin == 0
    load('tmp1.mat')
end


DATA = DefaultField(DATA,'REACTIONS_RESULTANTS_CALCULATE',[]) ;
DATA = DefaultField(DATA,'NODES_ENTITIES',[]) ;

ReactionResultantDOFS = [] ;
DATA.ReactionResultants = [] ;

if ~isempty(DATA.REACTIONS_RESULTANTS_CALCULATE)   && ~isempty(DATA.NODES_ENTITIES)
    DATA.REACTIONS_RESULTANTS_CALCULATE = DefaultField(DATA.REACTIONS_RESULTANTS_CALCULATE,'ENTITY',1) ;
    DATA.REACTIONS_RESULTANTS_CALCULATE = ...
        DefaultField(DATA.REACTIONS_RESULTANTS_CALCULATE,'DIRECTION','x') ;
    
    ReactionResultantDOFS = cell(1,length(DATA.REACTIONS_RESULTANTS_CALCULATE));
      DATA.OPERATORreactions = ReactionResultantDOFS ; % Not used so far (12-May-2020)
    for ientityLOC = 1:length(DATA.REACTIONS_RESULTANTS_CALCULATE.ENTITY)
        ientity = DATA.REACTIONS_RESULTANTS_CALCULATE.ENTITY(ientityLOC) ;
        
        
        NODES = DATA.NODES_ENTITIES{ientity} ;
        
        % Now we have to compute the degrees of freedom corresponding to
        % directions x and y 
        % --------------------------------------------------------------
        % 
        DOFS = cell2mat(DATA.TableDOFSnode(NODES)) ;   
         
        
        switch DATA.REACTIONS_RESULTANTS_CALCULATE.DIRECTION{ientityLOC}
            case 'x'
                DOFSloc = DOFS(1,:) ;
                
                
            case 'y'
                DOFSloc = DOFS(2,:) ;
                
            otherwise
                error('Option not implemented yet')
        end
        
        ReactionResultantDOFS{ientityLOC} = DOFSloc ;
        
    end
    
    DATA.ReactionResultantDOFS = ReactionResultantDOFS ;
    
    DATA.ReactionResultants = zeros(length(DATA.REACTIONS_RESULTANTS_CALCULATE.ENTITY),length(DATA.TIME_DISCRETIZATION));
    
    
end


