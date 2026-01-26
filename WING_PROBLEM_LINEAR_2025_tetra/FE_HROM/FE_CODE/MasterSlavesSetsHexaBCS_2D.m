function [MASTER SLAVES] =MasterSlavesSetsHexaBCS_2D(NODESln,NODESpnt,COOR,TOL,NormalPlanes)

% See HEXAG_2.jpg
if nargin == 0
    load('tmp1.mat')
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
 

%% POINTS
% ---------------
indNODMASTER.POINTS =   cell(1,1);
indNODSLAVES.POINTS =   cell(1,2);
NORMALS.POINTS =   cell(1,2) ;


imaster=0 ;
imaster = imaster + 1;
indNODMASTER.POINTS{imaster} = [] ;
indNODSLAVES.POINTS{imaster,1} = [1 ];
indNODSLAVES.POINTS{imaster,2} = [ 2 ];

% 
% imaster = imaster + 1;
% indNODMASTER.POINTS{imaster} = [2] ;
% indNODSLAVES.POINTS{imaster,1} = [4 ];
% indNODSLAVES.POINTS{imaster,2} = [ 6 ];
% indNODSLAVES.POINTS{imaster,3} = [ 8 ];
% indNODSLAVES.POINTS{imaster,4} = [ 10 ];
% indNODSLAVES.POINTS{imaster,5} = [ 12 ];

 NODESpnt = mat2cell(NODESpnt,ones(size(NODESpnt)),1); 

%%%%% SLAVE NODES


MASTER = cell(11,1) ;
SLAVES= cell(11,6) ;
FFF = fields(indNODMASTER) ;
imaster = 0 ;
for itype = 1:length(FFF)
    indNODMASTERloc = indNODMASTER.(FFF{itype}) ;
    indNODSLAVESloc = indNODSLAVES.(FFF{itype}) ;
    NORMALSloc = NORMALS.(FFF{itype}) ;
    switch FFF{itype}
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
        
        
        %
        %
        %         if length(indLOC)==2
        %             % Plane
        %             MASTER{imaster} = NODESpl{indLOC(1),indLOC(2)} ;
        %             SLAVES{imaster,1} = NODESpl{indLOC(1),indLOC(2)+1} ;
        %
        %             % Find exact correspondence between nodes
        %             % JAHOL --> Accommodate situations with more master nodes than
        %             % slaves nodes
        %             [SLAVES{imaster,1} MASTER{imaster}]= FindLocationMSlavePlane(MASTER{imaster},SLAVES{imaster,1},COOR,TOL,indLOC(1)) ;
        %             %
        %             %
        %
        %
        %         elseif length(indLOC)==4
        %             % Line
        %             MASTER{imaster} = NODESln{indLOC(1),indLOC(2),indLOC(3),indLOC(4)} ;
        %             SLAVES{imaster,1} = NODESln{indLOC(1),indLOC(2)+1,indLOC(3),indLOC(4)} ;
        %             SLAVES{imaster,2} = NODESln{indLOC(1),indLOC(2)+1,indLOC(3),indLOC(4)+1} ;
        %             SLAVES{imaster,3} = NODESln{indLOC(1),indLOC(2),indLOC(3),indLOC(4)+1} ;
        %
        %             for islave = 1:3
        %                 [ SLAVES{imaster,islave} MASTER{imaster}]= FindLocationMSlavePlane(MASTER{imaster},SLAVES{imaster,islave},COOR,TOL,[indLOC(1),indLOC(3)]) ;
        %             end
        %
        %         elseif length(indLOC)==1
        %             % Point
        %             %%% All vertices are regarded as MASTER NODES !!
        %             MASTER{imaster} = [] ;
        %             restP = NODESpnt ; %setdiff(NODESpnt,NODESpnt(indLOC(1))) ;
        %             for islave=1:length(restP)
        %                 SLAVES{imaster,islave}   = restP(islave) ;
        %             end
        %         end
        %     end
        
    end
end