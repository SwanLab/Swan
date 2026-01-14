function ReactionsPlot(DATA,OPERfe,NODES_SNAP)
%save('tmp1.mat')
if nargin == 0
    load('tmp1.mat')
    %     DATA.REACTIONS_RESULTANTS_CALCULATE.TIME_VARIABLE = 'FACTOR_TIME_DISPLACEMENTS';
    %     DATA.REACTIONS_RESULTANTS_CALCULATE.ENTITY = [3,4]                             ;
    %     DATA.REACTIONS_RESULTANTS_CALCULATE.DIRECTION = {'x','y'}                      ;
    %     DATA.REACTIONS_RESULTANTS_CALCULATE.TIME_VARIABLE = 'FACTOR_TIME_DISPLACEMENTS';
end

[~,NAMEFILE,~] = fileparts(DATA.nameWORKSPACE) ;

DATA = DefaultField(DATA,'REACTIONS_RESULTANTS_CALCULATE',[]) ;
DATA = DefaultField(DATA,'NODES_ENTITIES',[]) ;

DATA.REACTIONS_RESULTANTS_CALCULATE = DefaultField(DATA.REACTIONS_RESULTANTS_CALCULATE...
    ,'TIME_VARIABLE','TIME_DISCRETIZATION') ;

timeVAR = DATA.(DATA.REACTIONS_RESULTANTS_CALCULATE.TIME_VARIABLE) ;
COLORS = {'k','b','g','r','k*','m'} ;
xvar = {} ;
yvar = {} ;
COLORS_STORE = {} ;
hplot = [] ;
LEGG = {} ;
ifigure = 2700 ;
if ~isempty(DATA.REACTIONS_RESULTANTS_CALCULATE)   && ~isempty(DATA.NODES_ENTITIES)
    DATA.REACTIONS_RESULTANTS_CALCULATE = DefaultField(DATA.REACTIONS_RESULTANTS_CALCULATE,'ENTITY',1) ;
    DATA.REACTIONS_RESULTANTS_CALCULATE = DefaultField(DATA.REACTIONS_RESULTANTS_CALCULATE,'DIRECTION','norm') ;
    
    DATA = DefaultField(DATA,'ReactionResultants',[]) ;
    
    for ientityLOC = 1:length(DATA.REACTIONS_RESULTANTS_CALCULATE.ENTITY)
        ientity = DATA.REACTIONS_RESULTANTS_CALCULATE.ENTITY(ientityLOC) ;
        
        figure(ifigure) ;
        hold on
        grid on
        xlabel('Time Factor') ;
        ylabel('Force ') ;
        title(['Reaction forces '])
        LEGG{end+1} = [NAMEFILE,' Force ','(',DATA.REACTIONS_RESULTANTS_CALCULATE.DIRECTION{ientityLOC},')','  Entity =',...
            num2str(ientity)]  ;
        
        TIME_SHOW =timeVAR(DATA.STEPS_TO_STORE) ;
        
        
        if ~isempty(DATA.ReactionResultants)
            FORCE_TOTAL = DATA.ReactionResultants(ientityLOC,:) ;
            TIME_SHOW = timeVAR ;
        else
            
            NODES = DATA.NODES_ENTITIES{ientity} ;
            DOFS = small2large(NODES,OPERfe.ndim) ;
            
            
            switch DATA.REACTIONS_RESULTANTS_CALCULATE.DIRECTION{ientityLOC}
                case 'x'
                    DOFSloc = DOFS(1:OPERfe.ndim:end) ;
                    FORCE_TOTAL = sum(NODES_SNAP.Reactions(DOFSloc,:)) ;
                    
                case 'y'
                    DOFSloc = DOFS(2:OPERfe.ndim:end) ;
                    FORCE_TOTAL = sum(NODES_SNAP.Reactions(DOFSloc,:)) ;
               
                    error('Option not implemented yet')
            end
            
            
            
            
        end
        
        
        
        hplot(end+1 ) = plot(TIME_SHOW,FORCE_TOTAL,'Color',rand(1,3)) ;
        
        
        xvar{end+1} = TIME_SHOW ;
        yvar{end+1} = FORCE_TOTAL ;
        COLORS_STORE{end+1} = COLORS{ientityLOC} ;
        
        legend(hplot,LEGG)
        
    end
    
    
end


ReactionsPlot.xvar =xvar ;
ReactionsPlot.yvar =yvar ;
ReactionsPlot.COLORS_STORE =COLORS_STORE ;
ReactionsPlot.hplot =hplot ;
ReactionsPlot.LEGG =LEGG ;
ReactionsPlot.ifigure =ifigure ;

save(DATA.nameWORKSPACE,'ReactionsPlot','-append') ;
