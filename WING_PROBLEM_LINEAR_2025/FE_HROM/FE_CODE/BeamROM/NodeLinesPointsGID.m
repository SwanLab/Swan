function [NODES_POINTS,NODES_LINES] = NodeLinesPointsGID(SLICE_NAME)
if nargin == 0
    load('tmp2.mat')
end
% Read data from GID's file -1.dat (...Export calculation file),
% 4-July-2018, JAHO

nameDAT1 = [SLICE_NAME,'-1.dat'] ;
nameDAT2 = [SLICE_NAME,'-2.dat'] ;

if ~exist(nameDAT1,'file')
    [aaa,bbb]  =   fileparts(SLICE_NAME) ;
    nameDAT1 = [SLICE_NAME,'.gid',filesep,bbb,'-1.dat'] ;
    nameDAT2 = [SLICE_NAME,'.gid',filesep,bbb,'-2.dat'] ;
    
end
NODES_LINES = ExtractPropertiesGidDATFile(nameDAT1) ; 
NODES_POINTS = ExtractPropertiesGidDATFile(nameDAT2) ; 




%nameDAT2 = [SLICE_NAME,'-1.dat'] ;
% NODES_FACES = {} ;
% if exist(nameDAT1,'file')
%     DDD = load(nameDAT1) ;
%     if ~isempty(DDD)
%         LABELS_FACES = unique(DDD(:,2)) ;
%         nmaxLABEL = max(LABELS_FACES) ;
%         NODES_FACES = cell(1,nmaxLABEL) ;
%         for ifaceLOC = 1:length(LABELS_FACES)
%             iface = LABELS_FACES(ifaceLOC) ;
%             INDX = find(DDD(:,2) == iface) ;
%             NODES_FACES{iface} = DDD(INDX,1) ;
%         end
%     end
% end

end


function NODES_LINES = ExtractPropertiesGidDATFile(nameDAT1)

NODES_LINES = {} ;
if exist(nameDAT1,'file')
    DDD = load(nameDAT1) ;
    if ~isempty(DDD)
        LABELS_LINES = unique(DDD(:,2)) ;
        nmaxLABEL = max(LABELS_LINES) ;
        NODES_LINES = cell(1,nmaxLABEL) ;
        for ilineLOC = 1:length(LABELS_LINES)
            iline = LABELS_LINES(ilineLOC) ;
            INDX = find(DDD(:,2) == iline) ;
            NODES_LINES{iline} = DDD(INDX,1) ;
        end
    end
end
 
for i = 1:length(NODES_LINES)
    NODES_LINES{i} = unique(NODES_LINES{i}) ;
end

end