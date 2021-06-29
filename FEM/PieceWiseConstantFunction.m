classdef PieceWiseConstantFunction < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        integrator
        quadrature
    end
    
    properties (Access = private)
       mesh 
       fValues
       quadOrder
    end
    
    methods (Access = public)
        
        function obj = PieceWiseConstantFunction(cParams)
            obj.init(cParams);
        end
        
        function fNodal = projectToLinearNodalFunction(obj)
            obj.createQuadrature();
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
           LHS = obj.integrator.computeLHS();
        end
        
        function RHS = computeRHS(obj)
            fG = obj.computeFgauss;
            xG = obj.computeXgauss;
            RHS = obj.integrator.integrateFgauss(fG,xG,obj.quadOrder);                        
        end
        
    end
    
end