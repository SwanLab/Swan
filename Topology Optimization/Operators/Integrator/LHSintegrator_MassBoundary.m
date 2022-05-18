classdef LHSintegrator_MassBoundary < LHSintegrator

    properties (Access = private)
        quadType
    end

    methods (Access = public)
        
        function obj = LHSintegrator_MassBoundary(cParams)
            obj.init(cParams);
        end

        function LHS = compute(obj)
            LHS = obj.computeBoundaryMassMatrix();
        end
        
    end
    
    methods (Access = protected)
        
        function Mr = computeBoundaryMassMatrix(obj)
            s = obj.createIntegratorParams();
            nInt = numel(s.compositeParams);
            ndof = s.compositeParams{1}.dim.ndofs;
            LHS = sparse(ndof,ndof);
            for iInt = 1:nInt
                sL = s.compositeParams{iInt};
                sL.type     = 'MassMatrix';
                sL.quadType = 'LINEAR';
                lhs = LHSintegrator.create(sL);
                LHSadd = lhs.compute();
                LHS = LHS + LHSadd;
            end
            Mr = LHS;
        end
        
        function cParams = createIntegratorParams(obj)
            bMeshes  = obj.mesh.createBoundaryMesh();
            nBoxFaces = numel(bMeshes);
            d = obj.dim;
            for iMesh = 1:nBoxFaces
                bMesh = bMeshes{iMesh};
                m  = bMesh.mesh;
                s.dim  = d;
                s.mesh = m;
                s.globalConnec = bMesh.globalConnec;
                cParams.compositeParams{iMesh} = s;
            end
        end
        
    end
    
end