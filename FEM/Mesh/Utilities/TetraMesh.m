function m = TetraMesh(length, height, width, nx, ny, nz)
   cMeshGlobal =  HexaMesh(length, height, width, nx, ny, nz);
    Xiso     =  [-1 ,-1, -1;...
        1, -1, -1;...
        1, 1, -1;...
        -1, 1, -1;...
        -1, -1, 1;...
        1, -1, 1;...
        1, 1, 1;...
        -1, 1, 1;];
    connecIso    = delaunay(Xiso);
    ss.coord     = Xiso;
    ss.connec    = connecIso;
    localMesh    = Mesh.create(ss);
    nelem        = cMeshGlobal.nelem;
    bCutConnec   = cMeshGlobal.connec;
    connecIso    = localMesh.connec;
    nElemIso     = size(connecIso,1);
    nnodeSubMesh = size(connecIso,2);
    subConnec    = bCutConnec(:,connecIso');
    subConnec    = reshape(subConnec',[nnodeSubMesh,nelem*nElemIso])';
    sss.coord    = cMeshGlobal.coord;
    sss.connec   = subConnec;
    subMesh      = Mesh.create(sss);
    m = subMesh.computeCanonicalMesh();
end