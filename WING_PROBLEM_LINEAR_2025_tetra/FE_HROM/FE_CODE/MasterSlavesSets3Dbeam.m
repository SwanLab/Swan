function  [f1NOD f2NOD  NODESpl] =MasterSlavesSets3Dbeam(NODESpl,COOR,TOL,NormalPlanes,DATA) ;


% See HEXAG_2.jpg
if nargin == 0
    load('tmp2.mat')
end


%% PLANES (4 master, 4 slaves)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
indNODMASTER.PLANES =   cell(3,1);
indNODSLAVES.PLANES =   cell(3,1) ;
NORMALS.PLANES = cell(3,1) ;
imaster=0 ;
imaster = imaster + 1;
indNODMASTER.PLANES{imaster} = 1 ;
indNODSLAVES.PLANES{imaster} = 3 ;
NORMALS.PLANES{imaster} = NormalPlanes{1} ;


imaster = imaster + 1;
indNODMASTER.PLANES{imaster} = 2 ;
indNODSLAVES.PLANES{imaster} = 4 ;
NORMALS.PLANES{imaster} = NormalPlanes{2} ;

imaster = imaster + 1;
indNODMASTER.PLANES{imaster} = 5 ;
indNODSLAVES.PLANES{imaster} = 6 ;
NORMALS.PLANES{imaster} = NormalPlanes{5} ;


MASTER = cell(3,1) ;
SLAVES= cell(3,1) ;
FFF = fields(indNODMASTER) ;
imaster = 0 ;
for itype = 1:length(FFF)
    indNODMASTERloc = indNODMASTER.(FFF{itype}) ;
    indNODSLAVESloc = indNODSLAVES.(FFF{itype}) ;
    NORMALSloc = NORMALS.(FFF{itype}) ;
    %  switch FFF{itype}
    %     case 'LINES'
    NODES =  NODESpl;
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
%dbstop('70')
DATA.MakeMeshByRepetition = DefaultField(DATA.MakeMeshByRepetition,'DIRECTIONloc','z');
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
    case 'z'
        % dbstop('82')
         f2NOD= MASTER{3} ;
        f1NOD = SLAVES{3} ;
end

% We renumber NODESpl
if size(COOR,2) ==3
    NODESpl{1} =  MASTER{1} ; 
    NODESpl{3} =  SLAVES{1} ; 
    NODESpl{2} =  MASTER{2} ; 
    NODESpl{4} =  SLAVES{2} ; 
    NODESpl{5} =  MASTER{3} ; 
    NODESpl{6} =  SLAVES{3} ; 
end

