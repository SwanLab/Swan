function m = collapseMeshes(m1,m2)
s.coord    = [m1.coord;m2.coord];
s.connec   = [m1.connec;m2.connec+max(m1.connec(:))];
m          = Mesh.create(s);
[~,v,w]    = unique(m.coord,'rows','stable');
cc         = m.connec(:);
[~,~,w2]   = unique(cc);
connec2    = m.connec;
connec2(:) = v(w(w2));
ss.connec  = connec2;
ss.coord   = s.coord;
mm         = Mesh.create(ss);
m          = mm.computeCanonicalMesh();
end