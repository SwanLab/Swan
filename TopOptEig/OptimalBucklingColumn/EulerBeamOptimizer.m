classdef EulerBeamOptimizer < handle

    properties (Access = public)
        columnMesh
        value
        initValue
    end
    
    properties (Access = protected)
        optimizerType
    end
    
    properties (Access = private)
        cost
        mesh
        initValueType
        meshType
        designVariable
        eigenModes
        constraint
        optimizer
        nElem
        nConstraints
        columnLength
        nValues
        youngModulus
        inertiaMoment
        minThick
        maxThick
        maxIter
    end


    properties (Access = private)
        nIter
        mmaVarComputer
    end

    methods (Access = public)

        function obj = EulerBeamOptimizer(cParams)
            obj.init(cParams)
            obj.createMesh();
            obj.createDesignVariable();
            obj.createEigModes();
            obj.createCost();
            obj.createConstraint();
            obj.createOptimizer();
            obj.optimizer.solveProblem();
            obj.value = obj.designVariable.value;
            obj.PostProcess();
        end

    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.initValue     = cParams.value;
            obj.nElem         = 500;
            obj.nConstraints  = 3; 
            obj.columnLength  = 1; 
            obj.nValues       = obj.nElem+1;
            obj.youngModulus  = 1;
            obj.inertiaMoment = 1;  
            obj.optimizerType = 'fmincon'; %NullSpace';%'MMA';'AlternatingPrimalDual';%'fmincon'; % IPOPT';
            obj.initValueType = 'Random'; % Random/Constant/External Value
            obj.meshType      = 'Structured'; %Structured/Unstructured
            obj.maxIter       = 1000;
            obj.minThick(1:obj.nElem,1) = 0.5; %sqrt(0.5/pi);
            obj.minThick(obj.nElem+1)   = 0;
            obj.maxThick(1:obj.nElem,1) = 10; %sqrt(100/pi);
            obj.maxThick(obj.nElem+1)   = 10000;
        end
        
        function createMesh(obj)
            s.nElem = obj.nElem;
            s.columnLength = obj.columnLength;
            s.meshType = obj.meshType;
            s.type = 'LINE';
            m = MeshGenerator(s);
            obj.mesh = m.mesh;
        end

        function createDesignVariable(obj)
            s.initValue = obj.initValue;
            s.initValueType = obj.initValueType;
            s.type  = 'AreaColumn';
            s.mesh  = obj.mesh;
            des = DesignVariable.create(s);
            obj.designVariable = des;  
        end

        function createCost(obj)
            sF.type = 'firstEignValue_functional';            
            sC.weights = 1;
            sC.nShapeFuncs = 1;
            sC.designVar = obj.designVariable;         
            sC.shapeFuncSettings{1} = sF;
            obj.cost = Cost(sC);
        end
        
        function createEigModes(obj)
            s.mesh           = obj.mesh;
            s.inertiaMoment  = obj.inertiaMoment;
            s.youngModulus   = obj.youngModulus;
            s.designVariable = obj.designVariable;
            eigModes = EigModes(s);
            obj.eigenModes = eigModes;
        end

        function createConstraint(obj)
            sF1.eigModes       = obj.eigenModes;
            sF1.eigNum         = 1;
            sF1.type = 'doubleEig';  

            sF2.eigModes       = obj.eigenModes;  
            sF2.eigNum         = 2;
            sF2.type = 'doubleEig';  

            sF3.type = 'volumeColumn';    
            sF3.mesh = obj.mesh;

            sC.nShapeFuncs = 3;
            sC.designVar = obj.designVariable;   
            sC.dualVariable  = [];
            sC.shapeFuncSettings{1} = sF1;
            sC.shapeFuncSettings{2} = sF2;
            sC.shapeFuncSettings{3} = sF3;
            obj.constraint = Constraint(sC);
        end       

        function createOptimizer(obj)
            s = SettingsOptimizer();
            s.optimizerNames.type = obj.optimizerType;

            s.optimizerNames.primal = 'PROJECTED GRADIENT';
            s.uncOptimizerSettings.scalarProductSettings = obj.designVariable.scalarProduct;
            s.uncOptimizerSettings.designVariable   = obj.designVariable;
            s.monitoringDockerSettings.mesh = obj.mesh;
            s.monitoringDockerSettings.optimizerNames = s.optimizerNames;
            s.monitoringDockerSettings.refreshInterval = 1;
            s.designVar         = obj.designVariable;
            s.targetParameters.optimality_tol  = 0.0005; %obj.incrementalScheme.targetParams;
            s.targetParameters.constr_tol = 0.0005;
            s.cost              = obj.cost;
            s.constraint        = obj.constraint;
            s.incrementalScheme.iStep  = 1; %obj.incrementalScheme;
            s.incrementalScheme.nSteps = 1;
            sD.nConstraints = 3;
            s.dualVariable     = DualVariable(sD);              
            s.uncOptimizerSettings.ub = obj.maxThick;
            s.uncOptimizerSettings.lb = obj.minThick;        
            s.outputFunction.type        = 'Topology';
            s.outputFunction.iterDisplay = 'none';
            s.type = obj.optimizerType;
            s.outputFunction.monitoring  = MonitoringManager(s);                  
            s.maxIter           = obj.maxIter;
            s.constraintCase = {'INEQUALITY','INEQUALITY','INEQUALITY'};
            %s.primalUpdater = 'PROJECTED GRADIENT';

            obj.optimizer = Optimizer.create(s); 
        end


        function PostProcess(obj)
            s.designVariable = obj.designVariable;
            s.mesh           = obj.mesh;
            s.scale          = 1;
            post = PostProcessColumn(s);
            post.plotColumn();
            obj.columnMesh = post.m;
        end
        
    end
end