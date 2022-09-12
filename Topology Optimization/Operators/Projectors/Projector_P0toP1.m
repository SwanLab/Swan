classdef Projector_P0toP1 < handle

    properties (Access = private)
        mesh
        connec
        nelem
        nnode
        npnod
        quadOrder
        quadrature
        field
    end
    
    methods (Access = public)

        function obj = Projector_P0toP1(cParams)
            obj.init(cParams);
            obj.createQuadrature();
            obj.createField();
        end

        function xProj = project(obj, x)
            % x should be a P0Function, not a matrix!
            LHS = obj.computeLHS();
            RHS = obj.computeRHS(x);
            xProj = (LHS\RHS);
        end

    end

    methods (Access = private)
        
        function init(obj, cParams)
            obj.mesh   = cParams.mesh;
            obj.connec = cParams.connec;
            obj.quadOrder = 'LINEAR';
        end
        
        function q = createQuadrature(obj)
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature(obj.quadOrder);
            obj.quadrature = q;
        end

        function createField(obj)
            s.mesh               = obj.mesh;
            s.ndimf              = 1; 
            s.interpolationOrder = 'LINEAR';
            s.quadratureOrder    = 'QUADRATIC';
            obj.field = Field(s);
        end
        
        function LHS = computeLHS(obj)
            s.type  = 'MassMatrix';
            s.mesh  = obj.mesh;
            s.field = obj.field;
            lhs = LHSintegrator.create(s);
            LHS = lhs.compute();
        end
        
        function RHS = computeRHS(obj, f)
            s.type      = 'ShapeFunction';
            s.mesh      = obj.mesh;
%             s.meshType  = obj.mesh.type;
%             s.fType     = 'Gauss';
%             s.fGauss    = obj.computeFgauss();
%             s.xGauss    = obj.computeXgauss();
            s.quadOrder = obj.quadOrder;
            s.npnod     = obj.mesh.nnodes;
            s.globalConnec = obj.mesh.connec;
            rhs = RHSintegrator.create(s);

            fG = obj.computeFgauss(f);
            xG = obj.computeXgauss();
            RHS = rhs.computeFromFgauss(fG,xG);
        end
        
        function x = computeXgauss(obj)
            xG = obj.quadrature.posgp;
            x = repmat(xG,[1,1,obj.mesh.nelem]);
        end
        
        function f = computeFgauss(obj, f)
            % to be generalized for n dimensions
            % should be in P0Function
            ngaus = obj.quadrature.ngaus;
            fV(1,:) = f;
            f = repmat(fV,[ngaus,1]);
        end
        
    end

end

