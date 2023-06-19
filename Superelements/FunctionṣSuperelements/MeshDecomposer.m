function [out] = MeshDecomposer(in)
    subdomains = in.subdomains;
    mesh = in.mesh;
    type = in.type;

    subMesh = [];
    subBoundMesh = [];
    interMesh = [];

    xBaricenter = mesh.computeBaricenter();
    xBaricenter = xBaricenter(1,:);

    lengthDomain = max(mesh.coord(:,1));
    lengthSub = lengthDomain/subdomains;

    initSubCoord = 0;
    s.coord = mesh.coord;
    for i=1:subdomains
        isSub = (xBaricenter > initSubCoord) & (xBaricenter <= initSubCoord + lengthSub);
        subConnec = mesh.connec(isSub,:);
        s.connec = subConnec;
        subdomainMesh = Mesh(s);
        subMesh = cat(1,subMesh,subdomainMesh.computeCanonicalMesh());

        boundary = subMesh(i).createBoundaryMesh();
        subBoundMesh = cat(1,subBoundMesh,[boundary{1}, boundary{2}]);
        initSubCoord = initSubCoord + lengthSub;

        if type == "Three"
            interMesh = cat(1,interMesh,subBoundMesh(i,1)); % Interface mesh is equal to its left subdomain.
        end
    end

    out.subMesh = subMesh;
    out.subBoundMesh = subBoundMesh;
    out.interMesh = interMesh;
end