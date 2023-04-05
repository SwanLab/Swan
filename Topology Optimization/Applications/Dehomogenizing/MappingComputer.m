classdef MappingComputer < handle


    properties (Access = private)
        meshDisc
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
            LHS = obj.computeLHS();
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
            obj.meshDisc           = obj.mesh.createDiscontinuousMesh();
        end
        
        function computeLHS(obj)
            K = obj.computeStiffnessMatrix();
            In = obj.interpolator;
            Kn = In'*K*In;
            obj.LHS = Kn;
        end

        function K = computeStiffnessMatrix(obj)
            s.mesh = obj.mesh;
            s.type = 'StiffnessMatrix';
            s.fun  = P1DiscontinuousFunction.create(obj.mesh,1);
            lhs    = LHSintegrator.create(s);
            K      = lhs.compute();
        end

        function RHS = computeRHS(obj,iDim)
            aI = obj.dilatedOrientation{iDim};
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature('QUADRATIC');
            fG = aI.evaluate(q.posgp);            
            s.fType     = 'Gauss';
            s.fGauss    = fG;
            s.xGauss    = q.posgp;
            s.mesh      = obj.mesh;
            s.type      = obj.mesh.type;
            s.quadOrder = q.order;
            s.npnod     = obj.meshDisc.nnodes*1;
            s.type      = 'ShapeDerivative';
            s.globalConnec = obj.meshDisc.connec;
            rhs  = RHSintegrator.create(s);
            rhsV = rhs.compute();
            RHS = rhsV;
        end

        function u = solveSystem(obj,LHS,RHS)
            In = obj.interpolator;            
            LHS = In'*LHS*In;
            RHS = In'*RHS;            
            a.type = 'DIRECT';
            s = Solver.create(a);
            u = s.solve(LHS,RHS);
            u = In*u;            
            u = u(1:end);            
        end

    end

end
