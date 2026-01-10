function [ELEMENTS_SURFACES,NODES_LINES] = DAT_file_2dskeletonGID(SLICE_NAME)
% Reading ProblemType information: SKELETON2D.gid

nameDAT1 = [SLICE_NAME,'.dat'] ;
nameDAT2 = [SLICE_NAME,'-1.dat'] ;

if ~exist(nameDAT1,'file')
    [aaa,bbb]  =   fileparts(SLICE_NAME) ;
    nameDAT1 = [SLICE_NAME,'.gid',filesep,bbb,'.dat'] ;
    nameDAT2 = [SLICE_NAME,'.gid',filesep,bbb,'-1.dat'] ;
    
end


ELEMENTS_SURFACES = ExtractPropertiesDATfile(nameDAT1) ; 
NODES_LINES = ExtractPropertiesDATfile(nameDAT2) ; 

 


 
end 

function  NODES_FACES = ExtractPropertiesDATfile(nameDAT1)

NODES_FACES = {} ;
if exist(nameDAT1,'file')
    DDD = load(nameDAT1) ;
    if ~isempty(DDD)
        LABELS_FACES = unique(DDD(:,2)) ;
        nmaxLABEL = max(LABELS_FACES) ;
        NODES_FACES = cell(1,nmaxLABEL) ;
        for ifaceLOC = 1:length(LABELS_FACES)
            iface = LABELS_FACES(ifaceLOC) ;
            INDX = find(DDD(:,2) == iface) ;
            NODES_FACES{iface} = unique(DDD(INDX,1)) ;
        end
    end
else
    error('You have not generated the .dat file')
end

end