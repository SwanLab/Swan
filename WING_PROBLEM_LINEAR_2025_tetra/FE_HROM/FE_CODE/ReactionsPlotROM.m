function ReactionsFINAL = ReactionsPlotROM(DATAIN,NODES_SNAP,MESH2D,DOFsKEEP,NODES_LINES_ROTATIONS)

if nargin == 0
    load('tmp.mat')
elseif nargin==4
    NODES_LINES_ROTATIONS = [] ; 
end

DATAIN = DefaultField(DATAIN,'REACTIONS_RESULTANTS_CALCULATE',[]) ;
DATAIN = DefaultField(DATAIN,'NODES_ENTITIES',[]) ;
DATAIN.REACTIONS_RESULTANTS_CALCULATE = DefaultField(DATAIN.REACTIONS_RESULTANTS_CALCULATE...
    ,'TIME_VARIABLE','TIME_DISCRETIZATION') ;

%   RESULTS  at final time step
% -------------
ReactionsFINAL = zeros(DATAIN.ndimSP ,length( MESH2D.NODES_LINES ) ) ;
for ientity = 1:length( MESH2D.NODES_LINES )
    NODES = MESH2D.NODES_LINES{ientity} ;
    R = NODES_SNAP.Reactions(:,end) ;
    ResultantReactions = zeros(DATAIN.ndimSP,1) ; 
    
    for inode = 1:length(NODES) ; 
        inodeGLO = NODES(inode) ; 
        DOFS = DATAIN.TableDOFSnode{inodeGLO} ; 
        REACTloc = R(DOFS(1:DATAIN.ndimSP)) ; 
        if ~isempty(NODES_LINES_ROTATIONS)
        if ~isempty(NODES_LINES_ROTATIONS{ientity})
            if ~isempty(NODES_LINES_ROTATIONS{ientity}{inode})
                REACTloc = NODES_LINES_ROTATIONS{ientity}{inode}*REACTloc ; 
             end
        end
        end
        ResultantReactions =ResultantReactions + REACTloc ;  
    end
    
    ReactionsFINAL(:,ientity) = ResultantReactions ; 
    
       
     
    
    
    %             switch DATAIN.REACTIONS_RESULTANTS_CALCULATE.DIRECTION{ientityLOC}
    %                 case 'x'
    %                     DOFSloc = DOFS(1:ndim:end) ;
    %                     FORCE_TOTAL = sum(rALL(DOFSloc,:)) ;
    %                     if IS_NO_TIME == 1
    %                         FORCE_TOTAL = [0 FORCE_TOTAL] ;
    %                     end
    
    
    
    
end


timeVAR = DATAIN.(DATAIN.REACTIONS_RESULTANTS_CALCULATE.TIME_VARIABLE) ;
if isempty(timeVAR)
    timeVAR  = [0 1];
    IS_NO_TIME = 1;
else
    IS_NO_TIME = 0 ;
end


xvar = {} ;
yvar = {} ;
COLORS_STORE = {} ;
hplot = [] ;
LEGG = {} ;
ifigure = 2700 ;


COLORS = {'m-*','r','g','k','b','k-*','m'} ;

if ~isempty(DATAIN.REACTIONS_RESULTANTS_CALCULATE)
    DATAIN.REACTIONS_RESULTANTS_CALCULATE = DefaultField(DATAIN.REACTIONS_RESULTANTS_CALCULATE,'ENTITY',[1]) ;
    DATAIN.REACTIONS_RESULTANTS_CALCULATE = DefaultField(DATAIN.REACTIONS_RESULTANTS_CALCULATE,'DIRECTION',{'x'}) ;
    for ientityLOC = 1:length(DATAIN.REACTIONS_RESULTANTS_CALCULATE.ENTITY)
        ientity = DATAIN.REACTIONS_RESULTANTS_CALCULATE.ENTITY(ientityLOC) ;
        %   hplot = [] ;
        %   LEGG = [] ;
        figure(ifigure)
        hold on
        grid on
        xlabel('Time Factor')
        ylabel('Force ')
        title(['Reaction forces '])
        
        NAMELOLL = strrep(DATAIN.NAME_LOAD_DATAloc,'_',' ');
        
        DATAIN = DefaultField(DATAIN,'npointsECM',[]) ;
        if ~isempty(DATAIN.npointsECM)
            legPOINTS = ['m^* = ',num2str(DATAIN.npointsECM)] ;
        else
            legPOINTS = [] ;
        end
        
        
        LEGG{end+1} = [NAMELOLL,'Force ','(',DATAIN.REACTIONS_RESULTANTS_CALCULATE.DIRECTION{ientityLOC},')','  Entity =',...
            num2str(ientity), ' (ROM, DOFS =',num2str(DATAIN.ndimINTF),';',legPOINTS,')']  ;
        
        DATAIN = DefaultField(DATAIN,'ReactionResultants',[]) ;
        if ~isempty(DATAIN.ReactionResultants)
            FORCE_TOTAL = DATAIN.ReactionResultants(ientityLOC,:) ;
            TIME_SHOW = timeVAR ;
        else
            
            
            
            NODES = MESH2D.NODES_LINES{ientity} ;
            
            if ~isempty(DOFsKEEP)
                ndim= max(DATAIN.ndimINTF) ;
                ndof = size(MESH2D.COOR,1)*ndim ;
                R = NODES_SNAP.Reactions  ;
                rALL = zeros(ndof,size(R,2)) ;
                rALL(DOFsKEEP,:) = R ;
            else
                ndim =  DATAIN.ndimINTF(1) ;
                R = NODES_SNAP.Reactions  ;
                rALL = R;
            end
            
            DOFS = small2large(NODES,ndim) ;
            
            switch DATAIN.REACTIONS_RESULTANTS_CALCULATE.DIRECTION{ientityLOC}
                case 'x'
                    DOFSloc = DOFS(1:ndim:end) ;
                    FORCE_TOTAL = sum(rALL(DOFSloc,:)) ;
                    if IS_NO_TIME == 1
                        FORCE_TOTAL = [0 FORCE_TOTAL] ;
                    end
                    
                    
                case 'y'
                    DOFSloc = DOFS(2:ndim:end) ;
                    FORCE_TOTAL = sum(rALL(DOFSloc,:)) ;
                    if IS_NO_TIME == 1
                        FORCE_TOTAL = [0 FORCE_TOTAL] ;
                    end
                    
                otherwise
                    error('Option not implemented yet')
            end
            
        end
        hplot(end+1 ) = plot(timeVAR,FORCE_TOTAL,COLORS{ientityLOC}) ;
        legend(hplot,LEGG) ; 
        
        xvar{end+1} = timeVAR  ;
        yvar{end+1} = FORCE_TOTAL ;
        COLORS_STORE{end+1} = COLORS{ientityLOC} ;
        
    end
    
    
    ReactionsPlot.xvar =xvar ;
    ReactionsPlot.yvar =yvar ;
    ReactionsPlot.COLORS_STORE =COLORS_STORE ;
    ReactionsPlot.hplot =hplot ;
    ReactionsPlot.LEGG =LEGG ;
    ReactionsPlot.ifigure =ifigure ;
    
    DATAIN = DefaultField(DATAIN,'NameFile_msh',DATAIN.NAME_WS_MODES{1})  ;
    
    
    [FOLDER_PRINT,~,~] = fileparts(DATAIN.NameFile_msh) ;
    FOLDER_PRINT  = fileparts(FOLDER_PRINT) ;
    FOLDER_PRINT  = fileparts(FOLDER_PRINT) ;
    
    
    FOLDER_PRINT = [FOLDER_PRINT,filesep,'GRAPHS'] ;
    if ~exist(FOLDER_PRINT,'dir')
        mkdir(FOLDER_PRINT)
    end
    
    DATAIN_nameGRAPHS = [FOLDER_PRINT,filesep,DATAIN.NAME_LOAD_DATAloc,'.mat'] ;
    
    DATAIN = DefaultField(DATAIN,'nameGRAPHS',DATAIN_nameGRAPHS) ;
    
    if exist(DATAIN.nameGRAPHS) ==0
        save(DATAIN.nameGRAPHS,'ReactionsPlot','ReactionsFINAL') ;
    else
        save(DATAIN.nameGRAPHS,'ReactionsPlot','ReactionsFINAL','-append') ;
    end
    
    
end