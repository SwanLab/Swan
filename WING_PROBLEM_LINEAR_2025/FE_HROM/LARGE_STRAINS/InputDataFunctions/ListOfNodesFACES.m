function rnod = ListOfNodesFACES(NameFileMesh,ILINE,ndim)


try
    
    [NODES_FACES,NODES_LINES] = NodesFacesLinesGID(NameFileMesh(1:end-4)) ;
    
catch
    NAME_DATA = [NameFileMesh(1:end-4),'.gid',filesep,NameFileMesh(1:end-4)] ;
    [NODES_FACES,NODES_LINES] = NodesFacesLinesGID(NAME_DATA) ;
end

if ndim ==2
    rnod = unique(cell2mat(NODES_LINES(ILINE))) ;
    
else
    rnod = unique(cell2mat(NODES_FACES(ILINE))) ;
end