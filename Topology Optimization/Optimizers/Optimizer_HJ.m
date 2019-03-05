classdef Optimizer_HJ < Optimizer_Unconstrained
    
    properties (GetAccess = public, SetAccess = private)
        optimality_tol
    end
    
    properties (Access = private)
        e2
        meanCellSize
        filter
    end
    
    methods (Access = public)
        
        function obj = Optimizer_HJ(settings,epsilon,meanCellSize)
            obj@Optimizer_Unconstrained(settings,epsilon);
            obj.e2 = settings.e2;
            obj.meanCellSize = meanCellSize;
            obj.max_constr_change = +Inf;
            obj.nconstr = settings.nconstr;
            
            obj.setupFilter(settings,epsilon);
        end
        
        function phi = computeX(obj,phi,gradient)
            V = -obj.filter.regularize(phi,gradient);
            
            dt = 0.5*obj.e2*obj.line_search.kappa*obj.meanCellSize/max(abs(V(:))) ;
            phi = obj.solvelvlset(phi,V,dt);
            
            obj.opt_cond = obj.line_search.kappa;
        end
        
    end
    
    methods (Access = private)
        
        function solvedPhi = solvelvlset(obj,phi,V,dt)
            for i = 1:obj.line_search.HJiter
                phi = phi - dt*V;
            end
            solvedPhi = phi;
        end
        
        function setupFilter(obj,s,e)
            if obj.settingsFilterIsNotPDE(s)
                s.filter = 'PDE';
                obj.displayChangingFilter(s)
            end
            obj.filter = FilterFactory.create(s.filter,s.optimizer);
            obj.filter.setupFromGiDFile(s.filename,s.ptype);
            obj.filter.setDomainType('BOUNDARY');
            obj.filter.preProcess();
            obj.filter.updateEpsilon(e);
        end
        
    end
    
    methods (Static,Access = private)
        
        function itIsNot = settingsFilterIsNotPDE(s)
            if ~strcmp(s.filter,'PDE')
                itIsNot = true;
            else
                itIsNot = false;
            end
        end
        
        function displayChangingFilter(settings)
            print = settings.printChangingFilter;
            if isempty(print)
                print = true;
            end
            if print
                disp('Filter P1 changed to PDE for HJ velocity regularization');
            end
        end
        
    end
    
    methods
        
        function optimality_tol = get.optimality_tol(obj)
            optimality_tol = (0.0175/1e-3)*obj.target_parameters.optimality_tol;
        end
        
    end
    
end
