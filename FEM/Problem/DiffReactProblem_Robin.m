classdef DiffReactProblem_Robin < DiffReactProblem
    
    methods (Access = public)
        
        function computeVariables(obj,x)
            bc   = obj.boundaryConditions;
            xRed = bc.fullToReducedVector(x);
            LHS  = obj.computeLHS();
            x = obj.solver.solve(LHS,xRed);
            obj.variables.x = bc.reducedToFullVector(x);
        end
        
        function LHS = computeLHS(obj)
            Mr = obj.computeBoundaryMassMatrix();
            LHS = obj.epsilon^2*obj.K + obj.M + (obj.epsilon)*Mr;
            LHS = obj.boundaryConditions.fullToReducedMatrix(LHS);            
        end
        
    end

    methods (Access = private)
        
        function Mr = computeBoundaryMassMatrix(obj)
            s = obj.createIntegratorParams();
            nInt = numel(s.compositeParams);
            ndof = s.compositeParams{1}.dim.ndof;
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
            cParams.type  = 'COMPOSITE';
            cParams.npnod = obj.mesh.npnod;
            bMeshes  = obj.mesh.createBoundaryMesh();
            nBoxFaces = numel(bMeshes);
            d = obj.dim;
            for iMesh = 1:nBoxFaces
                bMesh = bMeshes{iMesh};
                m  = bMesh.mesh;
                s.dim  = d;
                s.mesh = m;
                s.type = 'SIMPLE';
                s.globalConnec = bMesh.globalConnec;
                cParams.compositeParams{iMesh} = s;
            end
        end
   
    end
    
end
