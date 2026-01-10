function [MASTER SLAVES] =MasterSlavesSetsHexa(NODESpl,NODESln,NODESpnt,COOR,TOL,NormalPlanes)

%dbstop('4')
if nargin == 0
    load('tmp1.mat')
end

%% PLANES (4 master, 4 slaves)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
indNODMASTER.PLANES =   cell(4,1);
indNODSLAVES.PLANES =   cell(4,1) ;
NORMALS.PLANES = cell(4,1) ;
imaster=0 ;
imaster = imaster + 1;
indNODMASTER.PLANES{imaster} = 1 ;
indNODSLAVES.PLANES{imaster} = 4 ;
NORMALS.PLANES{imaster} = NormalPlanes{1} ;


imaster = imaster + 1;
indNODMASTER.PLANES{imaster} = 2 ;
indNODSLAVES.PLANES{imaster} = 5 ;
NORMALS.PLANES{imaster} = NormalPlanes{2} ;

imaster = imaster + 1;
indNODMASTER.PLANES{imaster} = 3 ;
indNODSLAVES.PLANES{imaster} = 6 ;
NORMALS.PLANES{imaster} = NormalPlanes{3} ;

imaster = imaster + 1;
indNODMASTER.PLANES{imaster} = 7 ;
indNODSLAVES.PLANES{imaster} = 8 ;
NORMALS.PLANES{imaster} = NormalPlanes{7} ;

%% LINES
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
indNODMASTER.LINES =   cell(5,1);
indNODSLAVES.LINES =   cell(5,3) ;
NORMALS.LINES =   cell(5,3) ;

% Lines parallel to plane xy plane
%----------------------------------
imaster=0 ;
imaster = imaster + 1;
indNODMASTER.LINES{imaster} = 1 ;

indNODSLAVES.LINES{imaster,1} = [4 ];
NORMALS.LINES{imaster,1} =  NormalPlanes{1} ;

indNODSLAVES.LINES{imaster,2} = [7 ];
NORMALS.LINES{imaster,2} =  NormalPlanes{7} ;

indNODSLAVES.LINES{imaster,3} = [10];
NORMALS.LINES{imaster,3} =  NormalPlanes{1} + NormalPlanes{7};



imaster = imaster + 1;
indNODMASTER.LINES{imaster} = 2 ;

indNODSLAVES.LINES{imaster,1} = [5 ];
NORMALS.LINES{imaster,1} =  NormalPlanes{2} ;

indNODSLAVES.LINES{imaster,2} = [8 ];
NORMALS.LINES{imaster,2} =  NormalPlanes{7} ;

indNODSLAVES.LINES{imaster,3} = [11];
NORMALS.LINES{imaster,3} =  NormalPlanes{2} + NormalPlanes{7};


imaster = imaster + 1;
indNODMASTER.LINES{imaster} = 3 ;
indNODSLAVES.LINES{imaster,1} = [6 ];
NORMALS.LINES{imaster,1} =  NormalPlanes{3} ;

indNODSLAVES.LINES{imaster,2} = [9 ];
NORMALS.LINES{imaster,2} =  NormalPlanes{7} ;

indNODSLAVES.LINES{imaster,3} = [12];
NORMALS.LINES{imaster,3} =  NormalPlanes{3} + NormalPlanes{7};




% Lines parallel to  z  axis
imaster = imaster + 1;
indNODMASTER.LINES{imaster} = 13 ;
indNODSLAVES.LINES{imaster,1} = [15 ];
NORMALS.LINES{imaster,1} =  NormalPlanes{1} ;

indNODSLAVES.LINES{imaster,2} = [17 ];
NORMALS.LINES{imaster,2} =  NormalPlanes{2} ;


imaster = imaster + 1;
indNODMASTER.LINES{imaster} = 14 ;
indNODSLAVES.LINES{imaster,1} = [16 ];
NORMALS.LINES{imaster,1} =  NormalPlanes{2} ;

indNODSLAVES.LINES{imaster,2} = [18 ];
NORMALS.LINES{imaster,2} =  NormalPlanes{3} ;


%% POINTS
% ---------------
indNODMASTER.POINTS =   cell(2,1);
indNODSLAVES.POINTS =   cell(2,6);
NORMALS.POINTS =   cell(2,6) ;


imaster=0 ;
imaster = imaster + 1;
indNODMASTER.POINTS{imaster} = [] ;
indNODSLAVES.POINTS{imaster,1} = [1 ];
indNODSLAVES.POINTS{imaster,2} = [ 3 ];
indNODSLAVES.POINTS{imaster,3} = [ 5 ];
indNODSLAVES.POINTS{imaster,4} = [ 7 ];
indNODSLAVES.POINTS{imaster,5} = [ 9 ];
indNODSLAVES.POINTS{imaster,6} = [ 11];


imaster = imaster + 1;
indNODMASTER.POINTS{imaster} = [2] ;
indNODSLAVES.POINTS{imaster,1} = [4 ];
indNODSLAVES.POINTS{imaster,2} = [ 6 ];
indNODSLAVES.POINTS{imaster,3} = [ 8 ];
indNODSLAVES.POINTS{imaster,4} = [ 10 ];
indNODSLAVES.POINTS{imaster,5} = [ 12 ];

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
        case 'PLANES'
            NODES = NODESpl ;
        case 'LINES'
            NODES =  NODESln;
        case 'POINTS'
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