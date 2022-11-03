classdef EulerBeamOptimizer < handle

    properties (Access = public)
        columnMesh
        cost
    end
    
    properties (Access = protected)
        optimizerType
    end
    
    properties (Access = private)
        data
        bProblem
        nElem
        nConstraints
        columnLength
        mesh
        nValues
        youngModulus
        inertiaMoment
        minThick
        maxThick
        maxIter
    end

    properties (Access = private)
        dim
    end

    properties (Access = private)
        designVariable
        freeNodes
        nIter
        mmaVarComputer
    end

    methods (Access = public)

        function obj = EulerBeamOptimizer()
            obj.init()
            obj.createData();
            obj.createBucklingProblem();
            obj.createMMA();
            obj.computeIterativeProcess()
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
            obj.minThick      = 0.25;
            obj.maxThick      = 10;
            obj.optimizerType = 'MMA';
            obj.maxIter       = 100;
        end
        
        function createData(obj)
            s.nElem        = obj.nElem;
            s.columnLength = obj.columnLength;
            Data = DataCreator(s);
            obj.mesh = Data.mesh;
            obj.dim  = Data.dim;
            obj.freeNodes = Data.freeNodes;
            obj.data = Data;
        end
        
        function createBucklingProblem(obj)
            s.data = obj.data;
            s.youngModulus = obj.youngModulus;
            s.inertiaMoment = obj.inertiaMoment;
            problem = BucklingProblemCreator(s);
            obj.designVariable = problem.designVariable;
            obj.bProblem = problem;
        end       

        function obj = createMMA(obj)
            s.nConstraints  = obj.nConstraints;
            s.minThick      = obj.minThick;
            s.maxThick      = obj.maxThick;
            s.nValues       = obj.nValues;
            s.designVariable = obj.designVariable;
            s.mesh           = obj.mesh;
            mmaVarComp = MMAVariablesComputer(s);
            obj.mmaVarComputer = mmaVarComp;
        end

        function obj = computeIterativeProcess(obj)
            s.data = obj.data;
            s.bProblem = obj.bProblem;
            s.mmaVarComputer = obj.mmaVarComputer;
            s.nConstraints   = obj.nConstraints;
            s.nValues        = obj.nValues;
            s.youngModulus   = obj.youngModulus;
            s.inertiaMoment  = obj.inertiaMoment;
            s.maxIter        = obj.maxIter;
            s.optimizerType  = obj.optimizerType;
            s.designVariable = obj.designVariable;
            solution = IterativeProcessComputer(s);
            solution.compute();
            obj.cost = obj.designVariable.value(obj.nElem+1);
        end

        function PostProcess(obj)
            s.designVariable = obj.designVariable;
            s.dim            = obj.dim;
            s.mesh           = obj.mesh;
            s.scale          = 1;
            post = PostProcessColumn(s);
            post.plotColumn();
            obj.columnMesh = post.m;
        end
        
    end
end