classdef TractionLoad2 < handle
    
    properties (Access = public)
        type
        values
        fun
        mesh
    end
    
    methods (Access = public)
        
        function obj = TractionLoad2(mesh, s, type)
            obj.mesh = mesh;
            obj.type = type;
            
            switch type
                case 'DIRAC'
                    pl_dofs = s.domain(mesh.coord);
                    vals = zeros(length(pl_dofs), 1);
                    vals(pl_dofs, 1) = s.value;
                    obj.values = reshape(vals', [], 1);
                    
                case 'FUNCTION'
                    dom = s.domain;
                    f = s.fun;
                    neuFun = AnalyticalFunction.create(dom, f.mesh);
                    obj.fun = f .* neuFun;
            end
        end
        
        function rhs = computeRHS(obj,test)
            switch obj.type
                case 'DIRAC'
                    rhs = obj.values;
                case 'FUNCTION'
                    f   = obj.fun;
                    rhs = IntegrateRHS(@(v) DP(f,v),test,obj.mesh,'Boundary');
            end
        end
        
    end
    
end