function  [f1NOD f2NOD] =MasterSlavesSets2Dbeam(NODESln,COOR,TOL,NormalPlanes,DATA) ;


% See HEXAG_2.jpg
if nargin == 0
    load('tmp2.mat')
end


%% LINES
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
indNODMASTER.LINES =   cell(2,1);
indNODSLAVES.LINES =   cell(2,1) ;
NORMALS.LINES =   cell(2,1) ;

%----------------------------------
imaster=0 ;
imaster = imaster + 1;
indNODMASTER.LINES{imaster} = 1 ;

indNODSLAVES.LINES{imaster,1} = [3 ];
NORMALS.LINES{imaster,1} =  NormalPlanes{1} ;

imaster = imaster + 1;
indNODMASTER.LINES{imaster} = 2 ;

indNODSLAVES.LINES{imaster,1} = [4 ];
NORMALS.LINES{imaster,1} =  NormalPlanes{2} ;


% %% POINTS
% % ---------------
% indNODMASTER.POINTS =   cell(1,1);
% indNODSLAVES.POINTS =   cell(1,2);
% NORMALS.POINTS =   cell(1,2) ;

%
% imaster=0 ;
% imaster = imaster + 1;
% indNODMASTER.POINTS{imaster} = [] ;
% indNODSLAVES.POINTS{imaster,1} = [1 ];
% indNODSLAVES.POINTS{imaster,2} = [ 2 ];

%
% imaster = imaster + 1;
% indNODMASTER.POINTS{imaster} = [2] ;
% indNODSLAVES.POINTS{imaster,1} = [4 ];
% indNODSLAVES.POINTS{imaster,2} = [ 6 ];
% indNODSLAVES.POINTS{imaster,3} = [ 8 ];
% indNODSLAVES.POINTS{imaster,4} = [ 10 ];
% indNODSLAVES.POINTS{imaster,5} = [ 12 ];

%  NODESpnt = mat2cell(NODESpnt,ones(size(NODESpnt)),1);

%%%%% SLAVE NODES


MASTER = cell(2,1) ;
SLAVES= cell(2,1) ;
FFF = fields(indNODMASTER) ;
imaster = 0 ;
for itype = 1:length(FFF)
    indNODMASTERloc = indNODMASTER.(FFF{itype}) ;
    indNODSLAVESloc = indNODSLAVES.(FFF{itype}) ;
    NORMALSloc = NORMALS.(FFF{itype}) ;
    %  switch FFF{itype}
    %     case 'LINES'
    NODES =  NODESln;
    %    case 'POINTS'
    %   dbstop('111')
    %       NODES = NODESpnt  ;
    % end
    for imasterLOC = 1:length(indNODMASTERloc)
        imaster = imaster + 1 ;
        if ~isempty(indNODMASTERloc{imasterLOC})
            MASTER{imaster} = NODES{indNODMASTERloc{imasterLOC}} ;
        end
        for islave = 1:size(indNODSLAVESloc,2)
            %   if ~isempty(indNODSLAVESloc{imasterLOC,islave})
            SLAVESloc = NODES{indNODSLAVESloc{imasterLOC,islave}} ;
            if length(SLAVESloc) >1
                [SLAVES{imaster,islave} ]= FindLocationMSlavePlaneHEXAG(MASTER{imaster} ,SLAVESloc,COOR,TOL,NORMALSloc{imasterLOC,islave}) ;
            else
                SLAVES{imaster,islave} = SLAVESloc ;
            end
            %  else
            %     break
            % end
        end
        
        
    end
end

DATA.MakeMeshByRepetition = DefaultField(DATA.MakeMeshByRepetition,'DIRECTIONloc','x');
DIR = DATA.MakeMeshByRepetition.DIRECTIONloc ;

if iscell(DIR)
    DIR = DIR{1} ;
end

switch DIR
    case 'x'
        f1NOD = MASTER{1} ;
        f2NOD = SLAVES{1} ;
    case 'y'
        f1NOD = MASTER{2} ;
        f2NOD = SLAVES{2} ;
end

% % Therefore,
% f1NOD = MASTER{1} ;
% f2NOD = SLAVES{1} ;