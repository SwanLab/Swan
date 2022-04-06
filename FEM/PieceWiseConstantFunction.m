classdef PieceWiseConstantFunction < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        integrator
        quadrature
        dim
    end
    
    properties (Access = private)
       mesh 
       fValues
       quadOrder
    end
    
    methods (Access = public)
        
        function obj = PieceWiseConstantFunction(cParams)
            obj.init(cParams);
            obj.createQuadrature();
            obj.createDimensions();
        end
        
        function fNodal = projectToLinearNodalFunction(obj)
            obj.createIntegrator();
            LHS = obj.computeLHS();
            RHS = obj.computeRHS();
            fNodal = (LHS\RHS);
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh      = cParams.mesh;
            obj.fValues   = cParams.fValues;
            obj.quadOrder = 'LINEAR';
        end
        
        function q = createQuadrature(obj)
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature(obj.quadOrder);
            obj.quadrature = q;
        end

        function createDimensions(obj)
            q = Quadrature();
            q = q.set(obj.mesh.type);
            s.mesh = obj.mesh;
            s.pdim = '1D';
            s.ngaus = q.ngaus;
            d = DimensionVariables(s);
            d.compute();
            obj.dim = d;
        end

        function createIntegrator(obj)
            s.type = 'SIMPLE';
            s.mesh = obj.mesh;
            s.npnod = obj.mesh.npnod;
            s.globalConnec = obj.mesh.connec;
            int = Integrator.create(s);
            obj.integrator = int;
        end
        
        function x = computeXgauss(obj)
            xG = obj.quadrature.posgp;
            x = repmat(xG,[1,1,obj.mesh.nelem]);
        end
        
        function f = computeFgauss(obj)
            ngaus = obj.quadrature.ngaus;
            fV(1,:) = obj.fValues;
            f = repmat(fV,[ngaus,1]);
        end
        
        function LHS = computeLHS(obj)
            s.mesh         = obj.mesh;
            s.globalConnec = obj.mesh.connec;
            s.npnod        = obj.mesh.npnod;
            s.type         = 'MassMatrix';
            s.dim          = obj.dim;
            s.quadType     = 'QUADRATIC';
            lhs = LHSintegrator.create(s);
            LHS = lhs.compute();
        end
        
        function RHS = computeRHS(obj)
            fG = obj.computeFgauss;
            xG = obj.computeXgauss;
            RHS = obj.integrator.integrateFgauss(fG,xG,obj.quadOrder);
        end
        
    end
    
end