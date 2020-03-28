function exampleForSubMesh


coord = [0 0; 1 0; 1 1; 0 1; 1 2; 0 2];
connec = [1 2 3 4; 4 3 5 6];
s.coord = coord;
s.connec = connec;

m = Mesh().create(s);
m.plot;


s.mesh = m;
s.lastNode = max(connec(:));


figure
subMesher = SubMesher();
subMesh = subMesher.computeSubMesh(s);
subMesh.plot();
end

