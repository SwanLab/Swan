classdef Optimizer_HJ < Optimizer_Unconstrained
    
    properties (GetAccess = public, SetAccess = private)
        optimality_tol
    end
    
    properties (GetAccess = public, SetAccess = protected)
        type = 'HAMILTON JACOBI'
    end
    
    properties (Access = private)
        filter
    end
    
    methods (Access = public)
        
        function obj = Optimizer_HJ(cParams)
            designVar = cParams.designVariable;
            
            obj@Optimizer_Unconstrained(cParams);
            obj.setupFilter(cParams.scalarProductSettings.epsilon,designVar);
        end
        
        function phi = compute(obj)
            phi = obj.designVariable.value;
            g   = obj.objectiveFunction.gradient;
            phi = obj.updateLevelSet(phi,g);
            obj.updateOptimalityConditionValue();
            obj.designVariable.update(phi);
        end
       
    end
    
    methods (Access = private)
        
        function updateOptimalityConditionValue(obj)
            obj.optimalityCond = obj.lineSearch.value;
        end
        
        function V = computeVelocity(obj,g)
            V     = -obj.filter.regularize(g);
            Vnorm = max(abs(V(:)));
            V     = V/Vnorm;
        end
        
        function solvedPhi = updateLevelSet(obj,phi,g)
            V     = obj.computeVelocity(g);
            dt    = obj.lineSearch.value;
          %  for i = 1:obj.lineSearch.HJiter
                phi = phi - dt*V;
                phi = phi/norm(phi);
         %   end
            solvedPhi = phi;
        end
        
        function setupFilter(obj,e,designVar)
            s = SettingsFilter('paramsFilter_PDE_Boundary.json');
            s.mesh = designVar.mesh;
            s.designVarType = designVar.type;
            s.quadratureOrder = 'LINEAR';
            s.femSettings.scale = 'MACRO';
            s.designVariable = designVar;
            obj.filter = Filter.create(s);
            obj.filter.updateEpsilon(e);
        end
        
    end
    
end
