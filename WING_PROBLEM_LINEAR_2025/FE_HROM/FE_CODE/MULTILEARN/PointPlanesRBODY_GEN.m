function [NODESbound,DATALIM] =  PointPlanesRBODY_GEN(COORabs,CNref,DATA)


%dbstop('5')
if nargin ==0
    load('tmp2.mat')
end

NODESbound = [] ;

ndim = size(COORabs,2)  ;
xmin = min(COORabs(:,1)) ; xmin = xmin(1) ;   % Boundaries of domain i (square)
xmax = max(COORabs(:,1)) ; xmax = xmax(1) ;
ymin = min(COORabs(:,2)) ; ymin = ymin(1) ;
ymax = max(COORabs(:,2)) ; ymax = ymax(1) ;
if ndim ==2
    zmin = [] ;   zmax = [] ;
    DATA = DefaultField(DATA,'TypeUnitCell','HEXAG_2D_SQUARE') ;
else
    zmin = min(COORabs(:,3)) ; zmin = zmin(1) ;
    zmax = max(COORabs(:,3)) ; zmax = zmax(1) ;
    DATA = DefaultField(DATA,'TypeUnitCell','HEXAHEDRA') ;
end
DATALIM.xmin = xmin ; 
DATALIM.xmax = xmax ; 
DATALIM.ymin = ymin ; 
DATALIM.ymax = ymax ; 
DATALIM.zmin = zmin ; 
DATALIM.zmax = zmax ; 

TOL = ChooseTolerance(CNref,COORabs) ;
DATA = DefaultField(DATA,'CalculateMasterSlaves',1) ;
switch  DATA.TypeUnitCell
    case 'HEXAG_2D_SQUARE'
        %     error('Option not implemented')
        [NODESln, NormalPlanes]=  DetermineePlanesPeriodicHexaBCS(COORabs,TOL,xmax,xmin,ymax,ymin,zmax,zmin,1:size(COORabs,1)) ;
        
        
        [NODESpnt LINESpnt ]= DeterminePointsPeriodicGEN_2D(NODESln) ;
        
        [MASTER SLAVES] =MasterSlavesSets_2D(NODESln,COORabs,TOL,NormalPlanes);
        
        %   NODESbound.PLANE = cell(4,1) ;
        NODESbound.PLANE = cell(4,1) ;
        
        NODESbound.LINES = cell(4,1) ;
        
        
        NODESbound.PLANE{1} = MASTER{1}(:) ;   % face 1, xmin
        NODESbound.PLANE{2} = MASTER{2}(:) ;   % face 2  , ymin
        NODESbound.PLANE{3} =  SLAVES{1}(:) ; % face 3, xmax
        NODESbound.PLANE{4} =  SLAVES{2}(:) ; % face 4,  ymax
        
        
        NODESbound.POINTS = [] ;
        %%%%%%%%%%%%%%%%%%
        %%
        for ipoints = 1:length(NODESpnt)
            NODESbound.LINES{ipoints} = NODESpnt{ipoints} ;
        end
    case 'HEXAHEDRA'
        [NODESpl, NormalPlanes]=  DetermineePlanesPeriodicHexaBCS(COORabs,TOL,xmax,xmin,ymax,ymin,zmax,zmin,1:size(COORabs,1)) ;
        
        [NODESln PLANESln ]= DetermineLinesPeriodicGEN(NODESpl) ;
        
        [NODESpnt LINESpnt ]= DeterminePointsPeriodicGEN(NODESln) ;
        
        %   NODESpl = RemoveLinesPlanesHEXA(NODESpl,NODESln,PLANESln) ;
        
        % Remove points
        %         for ipoint = 1:length(NODESpnt)
        %             for iline = 1:length(NODESln)
        %                 INDD = find(NODESln{iline} == NODESpnt{ipoint}) ;
        %                 NODESln{iline}(INDD) = [] ;
        %             end
        %         end
        [MASTER SLAVES] =MasterSlavesSets_HEXAHEDRA(NODESpl,NODESln,NODESpnt,COORabs,TOL,NormalPlanes);
        
        NODESbound.PLANE = cell(6,1) ;
        NODESbound.LINES = cell(12,1) ;
        
        NODESbound.POINTS = cell(8,1) ;
        
        
        NODESbound.PLANE{1} = MASTER{1}(:) ;   % face 1, xmin
        NODESbound.PLANE{2} = MASTER{2}(:) ;   % face 2  , ymin
        NODESbound.PLANE{3} =  SLAVES{1,1}(:) ; % face 3, xmax
        NODESbound.PLANE{4} =  SLAVES{2,1}(:) ; % face 4,  ymax
        NODESbound.PLANE{5} = MASTER{3}(:) ;   % face 5,   zmax
        NODESbound.PLANE{6} =   SLAVES{3,1}(:) ; % face 6, zmin
        %%%%%%%%%%%%%%%%%%
        NODESbound.LINES{1} =  MASTER{4}(:) ;   % Line 1
        NODESbound.LINES{2} =  MASTER{5}(:) ;   % Line 2
        NODESbound.LINES{3}=  SLAVES{4,1}(:);               % Line 3
        NODESbound.LINES{4} =  SLAVES{5,1}(:);             % Line 4
        NODESbound.LINES{5} =  SLAVES{4,2}(:);            % Line 5
        NODESbound.LINES{6} =    SLAVES{5,2}(:);            % Line  6
        NODESbound.LINES{7} = SLAVES{4,3}(:);              % Line 7
        NODESbound.LINES{8} =    SLAVES{5,3}(:);             % Line 8
        NODESbound.LINES{9} = MASTER{6}(:) ;    % Line 9
        NODESbound.LINES{10} =  SLAVES{6,1}(:);            % Line 10
        NODESbound.LINES{11} = SLAVES{6,2}(:);             % Line 11
        NODESbound.LINES{12}=  SLAVES{6,3}(:);           % line 12
        %%
        for ipoints = 1:length(NODESpnt)
            NODESbound.POINTS{ipoints} = NODESpnt{ipoints} ;
        end
        
        
    otherwise
        error('Option not implemented')
end