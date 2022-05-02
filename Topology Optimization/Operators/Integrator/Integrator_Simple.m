classdef Integrator_Simple < Integrator
    
    methods (Access = public)
        
        function obj = Integrator_Simple(cParams)
            obj.init(cParams)
            obj.globalConnec = cParams.globalConnec;
        end
       
        function rhs = integrate(obj,fNodal)
            s.type      = 'ShapeFunction';
            s.mesh      = obj.mesh;
            s.meshType  = obj.mesh.type;
            s.fType     = 'Nodal';
            s.fNodal    = fNodal;
            s.quadOrder = 'LINEAR';
            s.npnod     = obj.npnod;
            s.globalConnec = obj.globalConnec;
            RHS = RHSintegrator.create(s);
            rhs = RHS.compute();
        end
        
    end
    
end