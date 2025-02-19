classdef LHSIntegratorMassBoundary < handle

    properties (Access = private)
        mesh
    end

    methods (Access = public)
        
        function obj = LHSIntegratorMassBoundary(cParams)
            obj.mesh  = cParams.mesh;
        end

        function LHS = compute(obj)
            LHS = obj.computeBoundaryMassMatrix();
        end
        
    end
    
    methods (Access = protected)
        
        function Mr = computeBoundaryMassMatrix(obj)
            s = obj.createIntegratorParams();
            nInt = numel(s.compositeParams);
            ndof = length(s.compositeParams{1}.bMesh.nodesInBoxFaces);
            LHSg = sparse(ndof,ndof);
            for iInt = 1:nInt
                sL = s.compositeParams{iInt};
                m = sL.mesh;
                a.type = 'MassMatrix';
                if m.nelem ~= 0
                    a.mesh = m;


                    a.test  = LagrangianFunction.create(sL.mesh, 1, 'P1');
                    a.trial = LagrangianFunction.create(sL.mesh, 1, 'P1');
                    lhs = LHSIntegrator.create(a);
                    LHS = lhs.compute();

                    local2global(sL.mesh.connec(:)) = sL.bMesh.globalConnec(:);
                    [iLoc,jLoc,vals] = find(LHS); % !!! iLoc, jLoc should come from P1Fun
                    iGlob = local2global(iLoc);
                    jGlob = local2global(jLoc);

                    LHSadd = sparse(iGlob,jGlob,vals, ndof, ndof);
                    LHSg = LHSg + LHSadd;
                end
            end
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