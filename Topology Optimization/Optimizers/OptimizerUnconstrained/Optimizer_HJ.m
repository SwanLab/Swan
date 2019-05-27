classdef Optimizer_HJ < Optimizer_Unconstrained
    
    properties (GetAccess = public, SetAccess = private)
        optimality_tol
    end
    
    properties (GetAccess = public, SetAccess = protected)
        type = 'HAMILTON JACOBI'
    end
    
    properties (Access = private)
        e2
        meanCellSize
        filter
    end
    
    methods (Access = public)
        
        function obj = Optimizer_HJ(cParams)
            designVar = cParams.designVariable;
            
            obj@Optimizer_Unconstrained(cParams);
            obj.e2 = cParams.e2;
            obj.meanCellSize = designVar.mesh.computeMeanCellSize();
            
            obj.setupFilter(cParams.scalarProductSettings.epsilon,designVar);
        end
        
        function phi = compute(obj)
            phi      = obj.designVariable.value;
            gradient = obj.objectiveFunction.gradient;                        
            V = -obj.filter.regularize(phi,gradient);
            
            dt = 0.5*obj.e2*obj.lineSearch.kappa*obj.meanCellSize/max(abs(V(:))) ;
            phi = obj.solvelvlset(phi,V,dt);
            obj.designVariable.value = phi;
            obj.opt_cond = obj.lineSearch.kappa;
        end
        
    end
    
    methods (Access = private)
        
        function solvedPhi = solvelvlset(obj,phi,V,dt)
            for i = 1:obj.lineSearch.HJiter
                phi = phi - dt*V;
            end
            solvedPhi = phi;
        end
        
        function setupFilter(obj,e,designVar)
            s = SettingsFilter('paramsFilter_PDE_Boundary.json');
            s.designVar = designVar;
            s.quadratureOrder = 'LINEAR';   
            s.femSettings.scale = 'MACRO';
            obj.filter = FilterFactory().create(s);
            obj.filter.preProcess();
            obj.filter.updateEpsilon(e);
        end
        
    end
    
end
