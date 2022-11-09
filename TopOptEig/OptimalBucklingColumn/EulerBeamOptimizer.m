classdef EulerBeamOptimizer < handle

    properties (Access = public)
        columnMesh
        cost
    end
    
    properties (Access = protected)
        optimizerType
    end
    
    properties (Access = private)
        mesh
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

        function obj = EulerBeamOptimizer()
            obj.init()
            obj.createMesh();
            obj.createDesignVariable();
            obj.createEigModes();
            obj.createCost();
            obj.createConstraint();
            obj.createOptimizer();
            obj.optimizer.solveProblem();
            obj.cost = obj.designVariable.value(obj.nElem+1);
            obj.PostProcess();
        end

    end
    
    methods (Access = private)
        
        function init(obj)
            obj.nElem         = 500;
            obj.nConstraints  = 3; 
            obj.columnLength  = 1; 
            obj.nValues       = obj.nElem+1;
            obj.youngModulus  = 1;
            obj.inertiaMoment = 1;  
            obj.minThick      = 0.5;
            obj.maxThick      = 10;
            obj.optimizerType = 'NullSpace';%NullSpace';%'MMA'; %'MMA'; 'AlternatingPrimalDual';%'fmincon';'NullSpace'; % IPOPT';
            obj.maxIter       = 100;
        end
        
        function createMesh(obj)
            s.coord  = obj.createCoordinates();  
            s.connec = obj.createConnectivity();
            s.type = 'LINE';
            m = Mesh(s);
            obj.mesh = m;
        end

        function createDesignVariable(obj)
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

        function coord = createCoordinates(obj)
            nnod = obj.nElem + 1;
%              x = [0;rand(nnod-2,1);1]*obj.columnLength;
%              x = sort(x);
%               coord = x; 
             coord = linspace(0,obj.columnLength,nnod)';
        end

        function connec = createConnectivity(obj)
            nNode = 2;
            Tnod = zeros(obj.nElem,nNode);
            e = 1;
            for iElem = 1: obj.nElem
                Tnod(iElem,1) = e;
                e = e + 1;
                Tnod(iElem,2) = e;
            end
            connec = Tnod;
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
            s.constraintCase = {'INEQUALITY','INEQUALITY','EQUALITY'};
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