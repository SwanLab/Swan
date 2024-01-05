classdef TopOpt_Problem < handle
    
    properties (GetAccess = public, SetAccess = public)
        designVariable
        dualVariable
        cost
        constraint
        physicalProblem
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
            obj.createHomogenizedVarComputer(cParams);
            obj.createEquilibriumProblem(cParams);
            obj.createCostAndConstraint(cParams);
            obj.createDualVariable();
            obj.createOptimizer(cParams);
            obj.createVideoMaker(cParams);
        end

        function createEquilibriumProblem(obj,cParams)
            s = cParams.problemData.femData;
            obj.physicalProblem = FEM.create(s);
            obj.initPrincipalDirections();
        end

        function initPrincipalDirections(obj)
            if isempty(obj.designVariable.alpha)
                dim = obj.physicalProblem.getDimensions();
                nelem = obj.designVariable.mesh.nelem;
                ndim = dim.ndimf;
                alpha0 = zeros(ndim,nelem);
                alpha0(1,:) = 1;
                obj.designVariable.alpha = alpha0;
            end
        end

        function createOptimizer(obj,settings)
            epsilon = obj.incrementalScheme.targetParams.epsilon;
            obj.completeOptimizerSettings(settings);
            obj.computeBounds();
            obj.optimizerSettings.outputFunction.type        = 'Topology';
            obj.optimizerSettings.outputFunction.iterDisplay = 'none';
            obj.optimizerSettings.outputFunction.monitoring  = MonitoringManager(obj.optimizerSettings);
            obj.optimizerSettings.uncOptimizerSettings.scalarProductSettings.femSettings.epsilon = epsilon;
            obj.optimizer = Optimizer.create(obj.optimizerSettings);
        end

        function computeBounds(obj)
            switch obj.designVariable.type
                case 'Density'
                obj.optimizerSettings.ub = 1;
                obj.optimizerSettings.lb = 0;
                case 'Density&Bound'
                    obj.optimizerSettings.uncOptimizerSettings.ub = ones(length(obj.designVariable.value),1);
                    obj.optimizerSettings.uncOptimizerSettings.ub(end) = 1000;
                    obj.optimizerSettings.uncOptimizerSettings.lb = zeros(length(obj.designVariable.value),1);
                    obj.optimizerSettings.uncOptimizerSettings.lb(end) = -1000;
                otherwise

            end
        end

        function completeOptimizerSettings(obj,cParams)
            s = cParams.optimizerSettings;
            
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
%             obj.optimizer.saveMonitoring();
%             obj.optimizer.simulationPrinter.print();
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
            s.mesh = obj.mesh;
        %    s.targetParamsSettings.epsilonPerInitial = 10*s.targetParamsSettings.epsilonPerFinal;
            obj.incrementalScheme = IncrementalScheme(s);
        end
        
        function createDesignVariable(obj,cParams)
            s = cParams.designVarSettings;
            s.mesh = obj.mesh;
            s.scalarProductSettings.epsilon = obj.incrementalScheme.targetParams.epsilon;
            s.scalarProductSettings.mesh    = obj.mesh;

            % (19/12/2023): The future idea will be to destroy
            % LevelSerCreator and use GeometricalFunction
            sLs        = s.creatorSettings;
            sLs.ndim   = obj.mesh.ndim;
            sLs.coord  = obj.mesh.coord;
            sLs.type   = s.initialCase;
            lsCreator  = LevelSetCreator.create(sLs);
            phi        = lsCreator.getValue();
            switch s.type
                case 'Density'
                    value = 1 - heaviside(phi);
                case 'LevelSet'
                    value = phi;
            end
            ss.fValues = value;
            ss.mesh    = obj.mesh;
            s.fun      = P1Function(ss);

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
            cParams.costSettings.physicalProblem       = obj.physicalProblem;
            cParams.constraintSettings.physicalProblem = obj.physicalProblem;
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

            a.coord  = s.femData.mesh.coord;
            a.connec = s.femData.mesh.connec;
            innerMesh = Mesh(a);
            boundMesh = innerMesh.createBoundaryMesh();
            sM.backgroundMesh = innerMesh;
            sM.boundaryMesh   = boundMesh;
%             obj.mesh = UnfittedMesh(sM);
            obj.mesh = innerMesh;
%             sM.coord  = s.femData.mesh.coord;
%             sM.connec = s.femData.mesh.connec;
%             obj.mesh = Mesh_Total(sM);
%             obj.mesh.innerMeshOLD.setMasterSlaveNodes(s.femData.mesh.masterSlaveNodes);
        end
        
    end
    
end
