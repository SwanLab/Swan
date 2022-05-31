classdef LHSintegrator_MassBoundary < LHSintegrator

    properties (Access = private)
        quadType
        field
    end

    methods (Access = public)
        
        function obj = LHSintegrator_MassBoundary(cParams)
            obj.mesh  = cParams.mesh;
            obj.field = cParams.field;
        end

        function LHS = compute(obj)
            LHS = obj.computeBoundaryMassMatrix();
        end
        
    end
    
    methods (Access = protected)
        
        function Mr = computeBoundaryMassMatrix(obj)
            s = obj.createIntegratorParams();
            nInt = numel(s.compositeParams);
            ndof = s.compositeParams{1}.field.dim.ndofs;
            LHS = sparse(ndof,ndof);
            for iInt = 1:nInt
                sL = s.compositeParams{iInt};
                sL.type     = 'MassMatrix';
                lhs = LHSintegrator.create(sL);
                LHSadd = lhs.compute();
                LHS = LHS + LHSadd;
            end
            Mr = LHS;
        end
        
        function cParams = createIntegratorParams(obj)
            bMeshes  = obj.mesh.createBoundaryMesh();
            nBoxFaces = numel(bMeshes);
            for iMesh = 1:nBoxFaces
                bMesh = bMeshes{iMesh};
                m  = bMesh.mesh;
                s.mesh = m;
                s.globalConnec = bMesh.globalConnec;
                a.ndimf = obj.field.dim.ndimf;
                a.mesh = m;
                a.interpolationOrder = 'LINEAR';
                a.quadratureOrder = 'LINEAR';
                a.scale = 'MACRO';
                subField = Field(a);
                subField.dim.ndofs = obj.field.dim.ndofs;
                s.field = subField;
                s.field.connec = bMesh.globalConnec;
                cParams.compositeParams{iMesh} = s;
            end
        end
        
    end
    
end