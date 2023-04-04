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
        eigenModes
        constraint
        optimizer    
    end
    
    properties (Access = private)
        initValueType
        meshType
        desVarType
        designVariable
        sectionVariables
        nElem
        nConstraints
        columnLength
        maxVolume
        nValues
        youngModulus
        inertiaMoment
        minDesVar
        maxDesVar
        minCost
        maxCost
        maxIter
        makeGIF
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
            obj.createSectionVariables();
            obj.createEigModes();
            obj.createCost();
            obj.createConstraint();
            obj.createOptimizer();
            obj.optimizer.solveProblem();
            obj.value = obj.designVariable.value;
            %obj.postProcess();
        end

    end
    
    methods (Access = private)
        
        function init(obj)
            obj.nElem         = 500;
            obj.nConstraints  = 3; 
            obj.columnLength  = 20; 
            obj.maxVolume     = 20*pi; %(1e-2)*20^3;
            obj.nValues       = obj.nElem+1;
            obj.youngModulus  = 1;
            obj.inertiaMoment = 1;  
            obj.makeGIF       = 'N'; %'Y';'N'
            obj.optimizerType = 'fmincon'; %NullSpace';%'MMA';'AlternatingPrimalDual';%'fmincon'; % IPOPT';
            obj.initValueType = 'Constant'; % Random/Constant/External Value/Sinus
            obj.meshType      = 'Structured'; %Structured/Unstructured
            obj.desVarType    = 'RadiusColumn'; %AreaColumn/RadiusColumn/SquareColumn/RectangularColumn/HoleColumn/RectangularHoleColumn
            obj.maxIter       = 100;
            obj.minDesVar = 0; 
            obj.maxDesVar = 200000; 
            obj.minCost   = 0;
            obj.maxCost   = 10000;
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
            s.type  = obj.desVarType;
            s.mesh  = obj.mesh;
            des = DesignVariable.create(s);
            obj.designVariable = des;  
        end

        function createSectionVariables(obj)
            s.mesh = obj.mesh;
            s.designVariable = obj.designVariable;
            sectVar = SectionVariablesComputer.create(s);
            obj.sectionVariables = sectVar;
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
            s.mesh             = obj.mesh;
            s.inertiaMoment    = obj.inertiaMoment;
            s.youngModulus     = obj.youngModulus;
            s.sectionVariables = obj.sectionVariables;
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
            sF3.sectionVariables = obj.sectionVariables;
            sF3.maxVolume = obj.maxVolume;

            sC.nShapeFuncs = 3;
            sC.designVar = obj.designVariable;   
            sC.dualVariable  = [];
            sC.shapeFuncSettings{1} = sF1;
            sC.shapeFuncSettings{2} = sF2;
            sC.shapeFuncSettings{3} = sF3;
            obj.constraint = Constraint(sC);
        end       

        function createOptimizer(obj)
            s.optimizerType    = obj.optimizerType; 
            s.desVar           = obj.designVariable;  
            s.sectionVariables = obj.sectionVariables;
            s.mesh             = obj.mesh;
            s.cost             = obj.cost;
            s.constraint       = obj.constraint; 
            s.minDesVar        = obj.minDesVar;
            s.maxDesVar        = obj.maxDesVar;
            s.minCost          = obj.minCost;
            s.maxCost          = obj.maxCost;
            s.maxIter          = obj.maxIter;
            creator = OptimizerCreator(s);
            obj.optimizer = creator.optimizer;
        end


        function postProcess(obj)
            s.designVariable = obj.designVariable;
            s.sectionVariables = obj.sectionVariables;
            s.mesh           = obj.mesh;
            s.scale          = 1;
            s.optimizer      = obj.optimizer;
            s.makeGIF        = obj.makeGIF;
            post = PostProcessColumn(s);
            post.plotColumn();
            obj.columnMesh = post.m;
        end
        
    end
end