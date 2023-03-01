classdef OptimizerCreator < handle
    
    properties (Access = public)
        optimizer
    end
    
    properties (Access = private)
        optimizerType 
        designVariable  
        sectionVariables
        mesh
        cost  
        constraint 
        upperBound
        lowerBound 
        maxIter 
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = OptimizerCreator(cParams)
            obj.init(cParams)
            obj.create() 
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.optimizerType   = cParams.optimizerType;
            obj.designVariable  = cParams.desVar;
            obj.sectionVariables=cParams.sectionVariables;
            obj.mesh            = cParams.mesh;
            obj.cost            = cParams.cost;
            obj.constraint      = cParams.constraint;
            obj.upperBound      = cParams.ub;
            obj.lowerBound      = cParams.lb;
            obj.maxIter         = cParams.maxIter;
            
        end
        
        function create(obj)
            s = SettingsOptimizer();
            s.optimizerNames.type = obj.optimizerType;
            s.optimizerNames.primal = 'PROJECTED GRADIENT';
            s.uncOptimizerSettings.scalarProductSettings = obj.designVariable.scalarProduct;
            s.uncOptimizerSettings.designVariable   = obj.designVariable;
            s.monitoringDockerSettings.mesh = obj.mesh;
            s.monitoringDockerSettings.optimizerNames = s.optimizerNames;
            s.monitoringDockerSettings.refreshInterval = 1;
            s.designVar         = obj.designVariable;
            s.sectionVariables  = obj.sectionVariables;
            s.targetParameters.optimality_tol  = 0.0005; %obj.incrementalScheme.targetParams;
            s.targetParameters.constr_tol = 0.0005;
            s.cost              = obj.cost;
            s.constraint        = obj.constraint;
            s.incrementalScheme.iStep  = 1; %obj.incrementalScheme;
            s.incrementalScheme.nSteps = 1;
            sD.nConstraints = 3;
            s.dualVariable     = DualVariable(sD);              
            s.uncOptimizerSettings.ub = obj.upperBound;
            s.uncOptimizerSettings.lb = obj.lowerBound;        
            s.outputFunction.type        = 'Topology';
            s.outputFunction.iterDisplay = 'none';
            s.type = obj.optimizerType;
            s.outputFunction.monitoring  = MonitoringManager(s);                  
            s.maxIter           = obj.maxIter;
            s.constraintCase = {'INEQUALITY','INEQUALITY','INEQUALITY'};
            %s.primalUpdater = 'PROJECTED GRADIENT';
            obj.optimizer = Optimizer.create(s);
        end
        
    end
    
end