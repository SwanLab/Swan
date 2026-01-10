function [MASTER SLAVES] =MasterSlavesSets_2D(NODESpl,COOR,TOL,NormalPlanes)

if nargin == 0
    load('tmp5.mat')
end

%% PLANES (3 master, 3 slaves)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
indNODMASTER.PLANES =   cell(2,1);
indNODSLAVES.PLANES =   cell(2,1) ;
NORMALS.PLANES = cell(2,1) ;
imaster=0 ;
imaster = imaster + 1;
indNODMASTER.PLANES{imaster} = 1 ;
indNODSLAVES.PLANES{imaster} = 3 ;
NORMALS.PLANES{imaster} = NormalPlanes{1} ;


imaster = imaster + 1;
indNODMASTER.PLANES{imaster} = 2 ;
indNODSLAVES.PLANES{imaster} = 4 ;
NORMALS.PLANES{imaster} = NormalPlanes{2} ;

 




MASTER = cell(2,1) ;
SLAVES= cell(2,1) ;
FFF = fields(indNODMASTER) ;
imaster = 0 ;
for itype = 1:length(FFF)
    indNODMASTERloc = indNODMASTER.(FFF{itype}) ;
    indNODSLAVESloc = indNODSLAVES.(FFF{itype}) ;
    NORMALSloc = NORMALS.(FFF{itype}) ;
    switch FFF{itype}
        case 'PLANES'
            NODES = NODESpl ;
%         case 'LINES'
%             NODES =  NODESln;
%         case 'POINTS'
%             %   dbstop('111')
%             NODES = NODESpnt  ;
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