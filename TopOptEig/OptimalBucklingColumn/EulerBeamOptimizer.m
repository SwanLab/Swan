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
        sectionVariables
        eigenModes
        constraint
        optimizer
        nElem
        nConstraints
        columnLength
        maxVolume
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
            obj.createSectionVariables();
            obj.createEigModes();
            obj.createCost();
            obj.createConstraint();
            obj.createOptimizer();
            obj.optimizer.solveProblem();
            obj.value = obj.designVariable.value;
            obj.postProcess();
        end

    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.initValue     = cParams.value;
            obj.nElem         = 500;
            obj.nConstraints  = 3; 
            obj.columnLength  = 20; 
            obj.maxVolume     = 80; %(1e-2)*20^3;
            obj.nValues       = obj.nElem+1;
            obj.youngModulus  = 1;
            obj.inertiaMoment = 1;  
            obj.optimizerType = 'fmincon'; %NullSpace';%'MMA';'AlternatingPrimalDual';%'fmincon'; % IPOPT';
            obj.initValueType = 'Random'; % Random/Constant/External Value
            obj.meshType      = 'Structured'; %Structured/Unstructured
            obj.maxIter       = 100;
            obj.minThick(1:obj.nElem,1) = sqrt(0.00001/pi); %sqrt(0.5/pi)/0.5/sqrt(0.5);
            obj.minThick(obj.nElem+1)   = 0;
            obj.maxThick(1:obj.nElem,1) = sqrt(1000/pi); %sqrt(100/pi)/10/sqrt(10);
            obj.maxThick(obj.nElem+1)   = 10000;
%             obj.minThick(1:2*obj.nElem,1) = 0.0005; 
%             obj.minThick(2*obj.nElem+1,1)   = 0;
%             obj.maxThick(1:2*obj.nElem,1) = 15; 
%             obj.maxThick(2*obj.nElem+1,1)   = 10000;
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
            s.type  = 'RadiusColumn'; %AreaColumn/RadiusColumn/SquareColumn/RectangularColumn/HoleColumn/RectangularHoleColumn
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
              s.optimizerType   = obj.optimizerType; 
              s.desVar          = obj.designVariable;  
              s.sectionVariables=obj.sectionVariables;
              s.mesh            = obj.mesh;
              s.cost            = obj.cost;
              s.constraint      = obj.constraint; 
              s.ub              = obj.maxThick;
              s.lb              = obj.minThick;
              s.maxIter         = obj.maxIter;
              creator = OptimizerCreator(s);
              obj.optimizer = creator.optimizer;
        end


        function postProcess(obj)
            s.designVariable = obj.designVariable;
            s.sectionVariables = obj.sectionVariables;
            s.mesh           = obj.mesh;
            s.scale          = 1;
            s.optimizer      = obj.optimizer;
            post = PostProcessColumn(s);
            post.plotColumn();
            obj.columnMesh = post.m;
        end
        
    end
end