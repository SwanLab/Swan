classdef Integrator_Simple < Integrator
    
    methods (Access = public)
        
        function obj = Integrator_Simple(cParams)
            obj.init(cParams)
            obj.globalConnec = cParams.globalConnec;
        end
       
        function rhs = integrate(obj,fNodal)
            quadOrder = 'LINEAR';
            rhs = obj.integrateFnodal(fNodal,quadOrder);
        end
        
        function rhs = integrateFnodal(obj,fNodal,quadOrder)
            s.type      = 'ShapeFunction';
            s.mesh      = obj.mesh;
            s.meshType  = obj.mesh.type;
            s.fType     = 'Nodal';
            s.fNodal    = fNodal;
            s.quadOrder = quadOrder;
            s.npnod     = obj.npnod;
            s.globalConnec = obj.globalConnec;
            RHS = RHSintegrator.create(s);
            rhs = RHS.compute();
        end

        % integrateFgauss is apparently used by PieceWiseConstantFunction
        function rhs = integrateFgauss(obj,fGauss,xGauss,quadOrder)
            type     = obj.mesh.type;
            rhsCells = obj.computeElementalRHS(fGauss,xGauss,type,quadOrder);
            rhs = obj.assembleIntegrand(rhsCells);
        end
        
    end
    
    methods (Access = private)
        
        function xGauss = computeGaussPoints(obj,quadOrder)
            q = obj.computeQuadrature(quadOrder);
            xGauss = repmat(q.posgp,[1,1,obj.mesh.nelem]);
        end
        
    end
    
end