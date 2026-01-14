function [MASTER SLAVES] =MasterSlavesSets_HEXAHEDRA(NODESpl,NODESln,NODESpnt,COOR,TOL,NormalPlanes)

if nargin == 0
    load('tmp5.mat')
end

%% PLANES (3 master, 3 slaves)
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



%% LINES
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
indNODMASTER.LINES =   cell(3,1);
indNODSLAVES.LINES =   cell(3,3) ;
NORMALS.LINES =   cell(3,3) ;

% Lines parallel to plane xy plane
%----------------------------------
imaster=0 ;
imaster = imaster + 1;
indNODMASTER.LINES{imaster} = 1 ;

indNODSLAVES.LINES{imaster,1} = [3 ];
NORMALS.LINES{imaster,1} =  NormalPlanes{1} ;

indNODSLAVES.LINES{imaster,2} = [5 ];
NORMALS.LINES{imaster,2} =  NormalPlanes{5} ;

indNODSLAVES.LINES{imaster,3} = [7];
NORMALS.LINES{imaster,3} =  NormalPlanes{1} + NormalPlanes{5};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imaster = imaster + 1;
indNODMASTER.LINES{imaster} = 2 ;

indNODSLAVES.LINES{imaster,1} = [4 ];
NORMALS.LINES{imaster,1} =  NormalPlanes{2} ;

indNODSLAVES.LINES{imaster,2} = [6 ];
NORMALS.LINES{imaster,2} =  NormalPlanes{5} ;

indNODSLAVES.LINES{imaster,3} = [8];
NORMALS.LINES{imaster,3} =  NormalPlanes{2} + NormalPlanes{5};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imaster = imaster + 1;
indNODMASTER.LINES{imaster} = 9 ;

indNODSLAVES.LINES{imaster,1} = [10 ];
NORMALS.LINES{imaster,1} =  NormalPlanes{2} ;

indNODSLAVES.LINES{imaster,2} = [11 ];
NORMALS.LINES{imaster,2} =  NormalPlanes{1} ;

indNODSLAVES.LINES{imaster,3} = [12];
NORMALS.LINES{imaster,3} =  NormalPlanes{2} + NormalPlanes{1};



% %% POINTS
% % ---------------
% indNODMASTER.POINTS =   cell(1,1);
% indNODSLAVES.POINTS =   cell(1,3);
% NORMALS.POINTS =   cell(1,4) ;


% imaster=0 ;
% imaster = imaster + 1;
% indNODMASTER.POINTS{imaster} = [] ;
% indNODSLAVES.POINTS{imaster,1} = [1 ];
% indNODSLAVES.POINTS{imaster,2} = [ 2 ];
% indNODSLAVES.POINTS{imaster,3} = [3 ];
% indNODSLAVES.POINTS{imaster,4} = [ 4 ];

%
% imaster = imaster + 1;
% indNODMASTER.POINTS{imaster} = [2] ;
% indNODSLAVES.POINTS{imaster,1} = [4 ];
% indNODSLAVES.POINTS{imaster,2} = [ 6 ];
% indNODSLAVES.POINTS{imaster,3} = [ 8 ];
% indNODSLAVES.POINTS{imaster,4} = [ 10 ];
% indNODSLAVES.POINTS{imaster,5} = [ 12 ];

% NODESpnt = mat2cell(NODESpnt,ones(size(NODESpnt)),1);

%%%%% SLAVE NODES


MASTER = cell(6,1) ;
SLAVES= cell(6,3) ;
FFF = fields(indNODMASTER) ;
imaster = 0 ;
for itype = 1:length(FFF)
    indNODMASTERloc = indNODMASTER.(FFF{itype}) ;
    indNODSLAVESloc = indNODSLAVES.(FFF{itype}) ;
    NORMALSloc = NORMALS.(FFF{itype}) ;
    switch FFF{itype}
        case 'PLANES'
            NODES = NODESpl ;
        case 'LINES'
            NODES =  NODESln;
        case 'POINTS'
            %   dbstop('111')
            NODES = NODESpnt  ;
    end
    for imasterLOC = 1:length(indNODMASTERloc)
        imaster = imaster + 1 ;
        if ~isempty(indNODMASTERloc{imasterLOC})
            MASTER{imaster} = NODES{indNODMASTERloc{imasterLOC}} ;
        end
        for islave = 1:size(indNODSLAVESloc,2)
            if ~isempty(indNODSLAVESloc{imasterLOC,islave})
                SLAVESloc = NODES{indNODSLAVESloc{imasterLOC,islave}} ;
                if length(SLAVESloc) >1
                    [SLAVES{imaster,islave} ]= FindLocationMSlavePlaneHEXAG(MASTER{imaster} ,SLAVESloc,COOR,TOL,NORMALSloc{imasterLOC,islave}) ;
                else
                    SLAVES{imaster,islave} = SLAVESloc ;
                end
            else
                break
            end
        end
        
        
        
    end
end