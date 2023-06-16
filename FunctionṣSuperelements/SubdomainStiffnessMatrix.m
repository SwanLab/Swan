function [subK] = SubdomainStiffnessMatrix(subMesh,material)
    subdomains = length(subMesh);
    subK = cell(subdomains,1);
    for i=1:subdomains
        s.fValues = zeros(subMesh(i).nnodes, subMesh(i).ndim);
        s.mesh    = subMesh(i);
        p1 = P1Function(s);

        s.type     = 'ElasticStiffnessMatrix';
        s.mesh     = subMesh(i);
        s.fun      = p1;
        s.material = material;
        lhs = LHSintegrator.create(s);
        K = lhs.compute();
        subK(i) = {K};
    end
end