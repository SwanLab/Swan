function [subC] = SubdomainLagrangeMatrix(Cinputs)
    subMesh = Cinputs.subMesh;
    subBoundMesh = Cinputs.subBoundMesh;
    lambdaType = Cinputs.lambdaType;

    subdomains = size(subMesh,1);
    boundaries = size(subBoundMesh,2);
    subC = cell(subdomains,boundaries);
    for i=1:subdomains
        for j=1:boundaries
            ndim  = subBoundMesh(i,j).mesh.ndim;

            test  = P1Function.create(subBoundMesh(i,j).mesh,ndim); % disp
            if lambdaType == "P0"
                trial = P0Function.create(subBoundMesh(i,j).mesh,ndim); % lambda
            elseif lambdaType == "P1"
                trial = P1Function.create(subBoundMesh(i,j).mesh,ndim); % lambda
            end

            s.type    = 'MassMatrix';
            s.mesh    = subBoundMesh(i,j).mesh;
            s.test    = test;
            s.trial   = trial;
            lhs       = LHSintegrator.create(s);
            C         = lhs.compute();
            subC(i,j) = {C};
        end
    end
end
