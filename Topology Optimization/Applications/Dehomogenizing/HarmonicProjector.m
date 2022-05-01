classdef HarmonicProjector < handle

    properties (Access = private)
        dim
        massMatrix
        stiffnessMatrix
        reducedStiffnessMatrix
        LHS
        solver
    end

    properties (Access = private)
        mesh
        boundaryMesh
    end

    methods (Access = public)

        function obj = HarmonicProjector(cParams)
            obj.init(cParams);
            obj.createDimension();
            obj.computeMassMatrix();
            obj.computeStiffnessMatrix();
            obj.computeReducedStiffnessMatrix();
            obj.computeLHS();
            obj.createSolver()
        end

        function lambda0 = computeInitalLambda(obj)
            b    = obj.boundaryMesh;
            nInt = setdiff(1:obj.dim.nnodes,b);
            lambda0 = zeros(length(nInt),1);
        end

        function [v,lambda] = solveProblem(obj,vH)
            lhs    = obj.LHS;
            rhs    = obj.computeRHS(vH);
            sol    = obj.solver.solve(lhs,rhs);
            v      = sol(1:obj.dim.nnodes,1);
            lambda = sol(obj.dim.nnodes+1:end,1);
        end

        function optPrim = computePrimalOptimaility(obj,lambda,v,vH)
            M        = obj.massMatrix;
            Kred     = obj.reducedStiffnessMatrix;
            gradPrim = M*(v-vH) + (Kred'*lambda);
            optPrim  = norm(gradPrim);
        end

        function vH = projectByDual(obj,v)
            Kred = obj.reducedStiffnessMatrix;
            M    = obj.massMatrix;
            lhs = Kred*(M\Kred');
            rhs = Kred*v;
            lambda  = obj.solver.solve(lhs,rhs);
            vH = v - M\(Kred'*lambda);
        end

        function cost = computeCost(obj,v,vH)
            M    = obj.massMatrix;
            cost = (v-vH)'*M*(v-vH)/(v'*M*v);
        end

        function optDual = computeDualOptimality(obj,vH)
            Kred = obj.reducedStiffnessMatrix;
            optDual = norm(Kred*vH);
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh             = cParams.mesh;
            obj.boundaryMesh     = cParams.boundaryMesh;
        end

        function createDimension(obj)
            q = Quadrature();
            q = q.set(obj.mesh.type);
            s.mesh = obj.mesh;
            s.pdim = '1D';
            s.ngaus = q.ngaus;
            d = DimensionVariables(s);
            d.compute();
            obj.dim = d;
        end

        function computeMassMatrix(obj)
            s.mesh         = obj.mesh;
            s.globalConnec = obj.mesh.connec;
            s.type         = 'MassMatrix';
            s.dim          = obj.dim;
            s.quadType     = 'QUADRATIC';
            lhs = LHSintegrator.create(s);
            M = lhs.compute();
            %M = diag(sum(M));
            %M = eye(size(M,1));
            obj.massMatrix = M;
        end

        function K = computeStiffnessMatrix(obj)
            s.mesh         = obj.mesh;
            s.globalConnec = obj.mesh.connec;
            s.type         = 'StiffnessMatrix';
            s.dim          = obj.dim;
            lhs = LHSintegrator.create(s);
            K = lhs.compute();
            obj.stiffnessMatrix = K;
        end

        function computeReducedStiffnessMatrix(obj)
            b    = obj.boundaryMesh;
            nInt = setdiff(1:obj.dim.nnodes,b);
            K    = obj.stiffnessMatrix;
            Kred = K(nInt,:);
            obj.reducedStiffnessMatrix = Kred;
        end

        function createSolver(obj)
            s = Solver.create();
            obj.solver = s;
        end

        function  computeLHS(obj)
            Kred = obj.reducedStiffnessMatrix;
            M    = obj.massMatrix;
            Z    = obj.computeZeroFunction();
            lhs  = [M,Kred';Kred,Z];
            obj.LHS = lhs;
        end

        function Z = computeZeroFunction(obj)
            b    = obj.boundaryMesh;
            nInt = setdiff(1:obj.dim.nnodes,b);
            Z    = zeros(length(nInt),length(nInt));
        end

        function rhs = computeRHS(obj,fNod)
            % Untested but should work
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature('LINEAR');
            s.type      = 'ShapeFunction';
            s.mesh      = obj.mesh;
            s.meshType  = obj.mesh.type;
            s.fType     = 'Nodal';
            s.fNodal    = fNod;
            s.quadOrder = q.order;
            s.npnod     = obj.dim.nnodes;
            s.globalConnec = obj.globalConnec;
            RHS = RHSintegrator.create(s);
            rhs = RHS.compute();
            b = obj.boundaryMesh;
            nInt = setdiff(1:obj.dim.nnodes,b);
            Z   = zeros(length(nInt),1);
            rhs = [rhs;Z];
        end

    end

end