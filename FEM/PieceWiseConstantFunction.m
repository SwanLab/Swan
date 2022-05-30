classdef PieceWiseConstantFunction < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
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
        
        function LHS = computeLHS(obj)
            s.mesh         = obj.mesh;
            s.globalConnec = obj.mesh.connec;
            s.type         = 'MassMatrix';
            s.dim          = obj.dim;
            s.quadType     = 'QUADRATIC';
            lhs = LHSintegrator.create(s);
            LHS = lhs.compute();
        end
        
        function RHS = computeRHS(obj)

            fG = obj.computeFgauss();
            xG = obj.computeXgauss();
        %    RHS = obj.integrator.integrateFgauss(fG,xG,obj.quadOrder);


            % Untested but should NOT work
            s.type      = 'ShapeFunction';
            s.mesh      = obj.mesh;
            s.meshType  = obj.mesh.type;
            s.fType     = 'Gauss';
            s.fGauss    = obj.computeFgauss();
            s.xGauss    = obj.computeXgauss();
            s.quadOrder = obj.quadOrder;
            s.npnod     = obj.mesh.nnodes;
            s.globalConnec = obj.mesh.connec;
            rhs = RHSintegrator.create(s);
            RHS = rhs.compute();
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

        
    end
    
end
