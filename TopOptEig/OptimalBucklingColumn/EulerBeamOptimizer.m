classdef EulerBeamOptimizer < handle
    
    properties (Access = protected)
        optimizerType
    end
    
    properties (Access = private)
        nElem
        nConstraints
        length
        nValues
        youngModulus
        inertiaMoment
        minThick
        maxThick
        maxIter 
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
            obj.createDesignVariable();
            obj.createBoundaryConditions()
            obj.createMMA();
            obj.computeIterativeProcess()
        end

    end
    
    methods (Access = private)
        
        function init(obj)
            obj.nElem         = 10;
            obj.nConstraints  = 3; 
            obj.length        = 1/obj.nElem; 
            obj.nValues       = obj.nElem+1;
            obj.youngModulus  = 1;
            obj.inertiaMoment = 1;  
            obj.minThick      = 0.25;
            obj.maxThick      = 10;
            obj.optimizerType = 'MMA';
            obj.maxIter       = 1000;
        end

        function createDesignVariable(obj)
            N = obj.nElem;
            s.type  = 'AreaColumn';
            s.nElem = obj.nElem;
            des = DesignVariable.create(s);
            x0 = ones(N+1,1);               
            des.update(x0);
            obj.designVariable = des;  
        end

        function createBoundaryConditions(obj)
            N = obj.nElem;
            fixnodes = union([1,2], [2*N+1,2*N+2]);
            nodes = 1:2*N+2;
            free  = setdiff(nodes,fixnodes);
            obj.freeNodes = free;
        end        

        function obj = createMMA(obj)
            s.nElem         = obj.nElem;
            s.nConstraints  = obj.nConstraints; 
            s.length        = obj.length;
            s.youngModulus  = obj.youngModulus;
            s.inertiaMoment = obj.inertiaMoment;
            s.minThick      = obj.minThick;
            s.maxThick      = obj.maxThick;
            s.nValues       = obj.nValues;
            s.x0            = obj.designVariable.value;
            mmaVarComp = MMAVariablesComputer(s);
            obj.mmaVarComputer = mmaVarComp;
        end

        function obj = computeIterativeProcess(obj)
            s.mmaVarComputer = obj.mmaVarComputer;
            s.freeNodes      = obj.freeNodes;
            s.nConstraints   = obj.nConstraints; 
            s.length         = obj.length;
            s.nValues        = obj.nValues;
            s.youngModulus   = obj.youngModulus;
            s.inertiaMoment  = obj.inertiaMoment;
            s.minThick       = obj.minThick;
            s.maxThick       = obj.maxThick;
            s.maxIter        = obj.maxIter;
            s.optimizerType  = obj.optimizerType;
            s.designVariable = obj.designVariable;
            s.nElem          = obj.nElem;
            solution = IterativeProcessComputer(s);
            solution.compute();
            obj.designVariable = solution.designVariable;
        end     
        
    end
end