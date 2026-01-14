function ReactionsPlotLIN(DATA,DATAOUT)
%save('tmp1.mat')
if nargin == 0
    load('tmp1.mat')
%     DATA.REACTIONS_RESULTANTS_CALCULATE.TIME_VARIABLE = 'FACTOR_TIME_DISPLACEMENTS';
%     
%     DATA.REACTIONS_RESULTANTS_CALCULATE.ENTITY = [3,4];
%     DATA.REACTIONS_RESULTANTS_CALCULATE.DIRECTION = {'x','y'};
%     DATA.REACTIONS_RESULTANTS_CALCULATE.TIME_VARIABLE = 'FACTOR_TIME_DISPLACEMENTS';
    
    
end

[~,NAMEFILE,~] = fileparts(DATA.nameWORKSPACE) ;
ndim = size(DATAOUT.posgp,1) ,, 

DATA = DefaultField(DATA,'REACTIONS_RESULTANTS_CALCULATE',[]) ;
DATA = DefaultField(DATA,'NODES_ENTITIES',[]) ;

DATA.REACTIONS_RESULTANTS_CALCULATE = DefaultField(DATA.REACTIONS_RESULTANTS_CALCULATE...
    ,'TIME_VARIABLE','TIME_DISCRETIZATION') ;


timeVAR =[0 1] ;


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
        
        NODES = DATA.NODES_ENTITIES{ientity} ;
        DOFS = small2large(NODES, ndim) ;
        
        switch DATA.REACTIONS_RESULTANTS_CALCULATE.DIRECTION{ientityLOC}
            case 'x'
                DOFSloc = DOFS(1:ndim:end) ;
                FORCE_TOTAL = sum(DATAOUT.React(DOFSloc,:)) ;  
                FORCE_TOTAL = [0,FORCE_TOTAL] ; 
                hplot(end+1 ) = plot(timeVAR,FORCE_TOTAL,COLORS{ientityLOC}) ;                
            case 'y'
                DOFSloc = DOFS(2:ndim:end) ;
                FORCE_TOTAL = sum(DATAOUT.React(DOFSloc,:)) ;
                   FORCE_TOTAL = [0,FORCE_TOTAL] ; 
                hplot(end+1 ) = plot(timeVAR,FORCE_TOTAL,COLORS{ientityLOC})  ;                 
            otherwise
                error('Option not implemented yet')
        end
        
        
        xvar{end+1} = timeVAR ;
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
