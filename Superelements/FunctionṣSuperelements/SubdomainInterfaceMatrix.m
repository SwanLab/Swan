function [subD] = SubdomainInterfaceMatrix(Dinputs)
    subBoundMesh = Dinputs.subBoundMesh;
    interMesh = Dinputs.interMesh;
    lambdaType = Dinputs.lambdaType;

    interfaces = size(interMesh,1);
    interfaceBoundaries = 2;
    subD = cell(interfaces,interfaceBoundaries);
    for i=2:interfaces % There is no interface at Dirichlet boundary.
        for j=1:interfaceBoundaries 
            ndim  = subBoundMesh(i,j).mesh.ndim;
            numSubdomain = i+j-2;
            numSide = 3-j;

            if lambdaType == "P0"
                test  = P0Function.create(subBoundMesh(numSubdomain,numSide).mesh,ndim); 
            elseif lambdaType == "P1"
                test  = P1Function.create(subBoundMesh(numSubdomain,numSide).mesh,ndim); 
            end
            trial = P1Function.create(interMesh(i).mesh,ndim); 
            
            s.type    = 'MassMatrix';
            s.mesh    = interMesh(i).mesh;
            s.test    = test;
            s.trial   = trial;
            lhs       = LHSintegrator.create(s);
            D         = lhs.compute();
            subD(i,j) = {D};
        end
    end
end
