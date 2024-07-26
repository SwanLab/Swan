classdef MappingComputer < handle


    properties (Access = private)
        mesh
        interpolator
        dilatedOrientation
        field
    end

    methods (Access = public)

        function obj = MappingComputer(cParams)
            obj.init(cParams);
        end

        function uF = compute(obj)
            LHS = obj.computeStiffnessMatrix();
            for iDim = 1:obj.mesh.ndim
                RHS = obj.computeRHS(iDim);
                uC  = obj.solveSystem(LHS,RHS);
                uC  = reshape(uC,obj.mesh.nnodeElem,[]);
                uV(iDim,:,:) = uC;
            end
            s.mesh    = obj.mesh;
            s.fValues = uV;
            uF = P1DiscontinuousFunction(s);
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh               = cParams.mesh;
            obj.dilatedOrientation = cParams.dilatedOrientation;
            obj.interpolator       = cParams.interpolator;
        end

        function K = computeStiffnessMatrix(obj)
            s.mesh  = obj.mesh;
            s.type  = 'StiffnessMatrix';
            s.test  = P1DiscontinuousFunction.create(obj.mesh, 1);
            s.trial = P1DiscontinuousFunction.create(obj.mesh, 1);
            lhs = LHSintegrator.create(s);
            K = lhs.compute();
        end

        function RHS = computeRHS(obj,iDim)
            aI = obj.dilatedOrientation{iDim};
            aI = aI.project('P1D');
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature('QUADRATIC');
            s.mesh      = obj.mesh;
            s.quadratureOrder = q.order;
            s.type      = 'ShapeDerivative';
            test = P1DiscontinuousFunction.create(obj.mesh,1);
            rhs  = RHSintegrator.create(s);
            rhsV = rhs.compute(aI,test);
            In = obj.interpolator;
            RHS = In'*rhsV;          
        end

        function u = solveSystem(obj,LHS,RHS)
            In = obj.interpolator;            
            LHS = In'*LHS*In;
            a.type = 'DIRECT';
            s = Solver.create(a);
            u = s.solve(LHS,RHS);
            u = In*u;            
            u = u(1:end);            
        end

    end

end
