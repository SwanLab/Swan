classdef LHSintegrator_MassBoundary_albert < handle

    properties (Access = private)
        mesh
    end

    methods (Access = public)
        
        function obj = LHSintegrator_MassBoundary_albert(cParams)
            obj.mesh  = cParams.mesh;
        end

        function LHS = compute(obj,trial,test)
            LHS = obj.computeBoundaryMassMatrix(trial,test);
        end
        
    end
    
    methods (Access = protected)
        
        function Mr = computeBoundaryMassMatrix(obj,trial,test)
            LHSg = sparse(ndof,ndof);
                a.type = 'MassMatrix';
                    a.mesh = m;


                    a.test  = LagrangianFunction.create(sL.mesh, 1, 'P1');
                    a.trial = LagrangianFunction.create(sL.mesh, 1, 'P1');
                    lhs = LHSintegrator.create(a);
                    LHS = lhs.compute();

                    local2global(sL.mesh.connec(:)) = sL.bMesh.globalConnec(:);
                    [iLoc,jLoc,vals] = find(LHS); % !!! iLoc, jLoc should come from P1Fun
                    iGlob = local2global(iLoc);
                    jGlob = local2global(jLoc);

                    LHSadd = sparse(iGlob,jGlob,vals, ndof, ndof);
                    LHSg = LHSg + LHSadd;
            Mr = LHSg;
        end
        
        function cParams = createIntegratorParams(obj)
            bMeshes  = obj.mesh.createBoundaryMesh();
            nBoxFaces = numel(bMeshes);
            for iMesh = 1:nBoxFaces
                bMesh = bMeshes{iMesh};
                s.bMesh = bMesh;
                s.mesh  = bMesh.mesh;
                cParams.compositeParams{iMesh} = s;
            end
        end
        
    end
    
end