classdef TopOpt_Problem < handle
    
    properties (GetAccess = public, SetAccess = public)
        designVariable
        dualVariable
        cost
        constraint
        optimizer
        incrementalScheme
        optimizerSettings
    end
    
    properties (Access = private)
        mesh
        videoMaker
        homogenizedVarComputer
    end
    
    methods (Access = public)
        
        function obj = TopOpt_Problem(cParams)
            obj.createMesh(cParams);
            obj.createIncrementalScheme(cParams);
            obj.createDesignVariable(cParams);
            obj.createHomogenizedVarComputer(cParams)
            obj.createCostAndConstraint(cParams);
            obj.createDualVariable();
            obj.createOptimizer(cParams);
            obj.createVideoMaker(cParams);
        end
        
        function createOptimizer(obj,settings)
            obj.completeOptimizerSettings(settings);
            obj.optimizer = Optimizer.create(obj.optimizerSettings);
        end
        
        function completeOptimizerSettings(obj,cParams)
            s = cParams.optimizerSettings;
            s.uncOptimizerSettings.scalarProductSettings = obj.designVariable.scalarProduct;
            
            s.uncOptimizerSettings.targetParameters = obj.incrementalScheme.targetParams;
            s.uncOptimizerSettings.designVariable   = obj.designVariable;
            
            s.monitoringDockerSettings.mesh = obj.mesh;
            
            s.designVar         = obj.designVariable;
            s.targetParameters  = obj.incrementalScheme.targetParams;
            s.cost              = obj.cost;
            s.constraint        = obj.constraint;
            s.incrementalScheme = obj.incrementalScheme;
            s.dualVariable      = obj.dualVariable;
            
            obj.optimizerSettings = s;
        end
        
        function computeVariables(obj)
            while obj.incrementalScheme.hasNext()
                obj.incrementalScheme.next();
                obj.optimizer.solveProblem();
            end
            obj.optimizer.saveMonitoring();
        end
        
        function postProcess(obj)
            iter = 0:obj.optimizer.nIter;
            obj.videoMaker.iterations = iter;
            obj.videoMaker.makeVideo();
        end
        
    end
    
    methods (Access = private)
        
        function optSet = obtainOptimizersSettings(obj,settings)
            epsilon = obj.incrementalScheme.targetParams.epsilon;
            settings.optimizerSettings.uncOptimizerSettings.lineSearchSettings.epsilon = epsilon;
            settings.optimizerSettings.uncOptimizerSettings.scalarProductSettings.epsilon = epsilon;
            set = settings.clone();
            optSet = set.optimizerSettings;
        end
        
        function createIncrementalScheme(obj,cParams)
            s = cParams.incrementalSchemeSettings;
            s.mesh = obj.mesh.innerMeshOLD;
        %    s.targetParamsSettings.epsilonPerInitial = 10*s.targetParamsSettings.epsilonPerFinal;
            obj.incrementalScheme = IncrementalScheme(s);
        end
        
        function createDesignVariable(obj,cParams)
            s = cParams.designVarSettings;
            s.mesh = obj.mesh;
            s.scalarProductSettings.epsilon = obj.incrementalScheme.targetParams.epsilon;
            s.scalarProductSettings.mesh    = obj.mesh.innerMeshOLD;            
            obj.designVariable = DesignVariable.create(s);
        end
        
        function createDualVariable(obj)
            s.nConstraints = obj.constraint.nSF;
            obj.dualVariable = DualVariable(s);
        end
        
        function createHomogenizedVarComputer(obj,cParams)
            s = cParams.homogenizedVarComputerSettings;
            obj.homogenizedVarComputer = HomogenizedVarComputer.create(s);
        end
        
        function createCostAndConstraint(obj,cParams)
            obj.createCost(cParams);
            obj.createConstraint(cParams);
        end
        
        function createCost(obj,cParams)
            s = cParams.costSettings;
            s.designVar = obj.designVariable;
            s.homogenizedVarComputer = obj.homogenizedVarComputer;
            s.targetParameters = obj.incrementalScheme.targetParams;
            obj.cost = Cost(s);
        end
        
        function createConstraint(obj,cParams)
            s = cParams.constraintSettings;
            s.designVar = obj.designVariable;
            s.homogenizedVarComputer = obj.homogenizedVarComputer;
            s.targetParameters = obj.incrementalScheme.targetParams;
            s.dualVariable = obj.dualVariable;
            obj.constraint = Constraint(s);
        end
        
        function createVideoMaker(obj,cParams)
            s = cParams.videoMakerSettings;
            obj.videoMaker = VideoMaker.create(s);
        end
        
        function createMesh(obj,cParams)
            s = cParams.designVarSettings;
            sM.coord  = s.femData.mesh.coord;
            sM.connec = s.femData.mesh.connec;
            obj.mesh = Mesh_Total(sM);
            obj.mesh.innerMeshOLD.setMasterSlaveNodes(s.femData.mesh.masterSlaveNodes);
        end
        
    end
    
end
