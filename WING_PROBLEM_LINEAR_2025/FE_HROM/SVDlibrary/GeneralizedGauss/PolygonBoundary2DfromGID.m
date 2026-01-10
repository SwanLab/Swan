function COORpolygonBND = PolygonBoundary2DfromGID(DATALOC)



DATALOC = DefaultField(DATALOC,'NameFileMeshTRI',[]) ;
if ~isempty(DATALOC.NameFileMeshTRI)
    NameFileMesh  = DATALOC.NameFileMeshTRI;
    [COORloc,CNloc,~,CONNECTbound]=...
        ReadMeshFile(NameFileMesh,'READ_MATERIAL_COLUMN',1)  ;
    [AAAA,BBB] =fileparts(NameFileMesh);
    NAmeDAT = [NameFileMesh(1:end-4),'.gid',filesep,BBB] ;
    [NODES_FACES,NODES_LINES] = NodesFacesLinesGID(NAmeDAT) ;
    % Outer polygon
    POLYGONS_DOMAIN = [] ;
    iline = 1;
    OUTERnodes=  NODES_LINES{iline} ;   % List of nodes delimiting external boundary
    OUTERnodesSORTED = PoligonFromListOfNodes(OUTERnodes,CONNECTbound)     ;
    xPOLY = COORloc(OUTERnodesSORTED,1)' ;
    yPOLY = COORloc(OUTERnodesSORTED,2)' ;
    % INNER BOUNDARIES
    for ibnd = 2:length(NODES_LINES)
        xPOLY =[xPOLY, NaN] ;
        yPOLY =[yPOLY, NaN] ;
        INNERNODES = PoligonFromListOfNodes(NODES_LINES{ibnd},CONNECTbound)  ;
        xPOLY =[xPOLY,  COORloc(INNERNODES,1)' ] ;
        yPOLY =[yPOLY,  COORloc(INNERNODES,2)' ] ;
    end
    
%       figure
% % %     
%       plot(xPOLY,yPOLY,'LineWidth',2) % polygon
%       axis equal
%     
    COORpolygonBND = [xPOLY',yPOLY'];
    
    
else
    error('PROVIDE A COARSE gid MESH IN WHICH THE INNER AND OUTTER BOUNDARIES ARE IDENTIFIED. DATA_GENGAUSS.NameFileMeshTRI)')
    
end
