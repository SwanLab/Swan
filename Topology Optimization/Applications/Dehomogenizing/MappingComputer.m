classdef MappingComputer < handle


    properties (Access = private)
        meshDisc
        mesh
        orientationVector
        dilatedOrientation
        field
    end

    methods (Access = public)

        function obj = MappingComputer(cParams)
            obj.init(cParams);
            obj.createField();
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
            obj.orientationVector  = cParams.orientationVector;
            obj.dilatedOrientation = cParams.dilatedOrientation;                       
            obj.meshDisc           = obj.mesh.createDiscontinuousMesh();
        end

        function createField(obj)
            s.mesh               = obj.meshDisc;
            s.ndimf              = 1;
            s.interpolationOrder = obj.mesh.interpolation.order;
            obj.field = Field(s);
        end
        
        function K = computeLHS(obj)
            s.mesh         = obj.mesh;
            s.globalConnec = obj.mesh.connec;
            s.type         = 'StiffnessMatrix';
            s.field        = obj.field;
            lhs = LHSintegrator.create(s);
            K = lhs.compute();
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
            s.npnod     = obj.field.dim.ndofs;
            s.type      = 'ShapeDerivative';
            s.globalConnec = obj.meshDisc.connec;
            rhs  = RHSintegrator.create(s);
            rhsV = rhs.compute();
            RHS = rhsV;
        end

        function u = solveSystem(obj,LHS,RHS)
            In = obj.orientationVector.interpolator;            
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