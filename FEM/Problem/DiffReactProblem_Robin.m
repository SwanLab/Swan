classdef DiffReactProblem_Robin < DiffReactProblem
    
    methods (Access = public)
        
        function computeVariables(obj,x)
            LHS  = obj.computeLHS();
            xReg = obj.solver.solve(LHS,x);
            obj.variables.x = xReg;
        end
        
        function LHS = computeLHS(obj)
            Mr = obj.computeBoundaryMassMatrix();
            LHS = obj.epsilon^2*obj.K + obj.M + (obj.epsilon)*Mr;
        end
        
    end

    methods (Access = private)
        
        function Mr = computeBoundaryMassMatrix(obj)
            cParams = obj.createIntegratorParams();
            nInt = numel(cParams.compositeParams);
            ndof = cParams.compositeParams{1}.dim.ndof;
            LHS = sparse(ndof,ndof);
            for iInt = 1:nInt
                s = cParams.compositeParams{iInt};
                s.type     = 'MassMatrix';
                s.quadType = 'LINEAR';
                lhs = LHSintegrator.create(s);
                LHSadd = lhs.compute();
                LHS = LHS + LHSadd;
            end
            Mr = LHS;
        end
        
        function params = createIntegratorParams(obj)
            params.type  = 'COMPOSITE';
            params.npnod = obj.mesh.npnod;
            fileName = obj.problemData.fileName;
            bMeshes  = obj.createBoundaryMesh(fileName);
            nBoxFaces = numel(bMeshes);
            d = obj.dim;
            for iMesh = 1:nBoxFaces
                boxFaceMesh = bMeshes{iMesh};
                bfMesh  = boxFaceMesh.mesh;
                gConnec = boxFaceMesh.globalConnec;
                cParams.dim  = d;
                cParams.mesh = bfMesh;
                cParams.type = 'SIMPLE';
                cParams.globalConnec = gConnec;
                params.compositeParams{iMesh} = cParams;
            end
        end
        
        function boundaryMesh = createBoundaryMesh(obj,fileName)
            run(fileName);
            if exist('External_border_nodes','var') && ~isempty(External_border_nodes)
                s.borderNodes    = External_border_nodes;
                s.borderElements = External_border_elements;
                s.backgroundMesh = obj.mesh;
                s.type = 'FromData';
                b = BoundaryMeshCreator.create(s);
                boundaryMesh = b.create();
            else
                s.backgroundMesh = obj.mesh;
                s.dimension = 1:obj.mesh.ndim;
                s.type = 'FromReactangularBox';
                bC = BoundaryMeshCreator.create(s);
                boundaryMesh = bC.create();
            end
        end

    end
    
end
