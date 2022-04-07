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
        nIter
        mmaVarComputer
    end
     
    methods (Access = public)
        
        function obj = EulerBeamOptimizer()
            obj.init()
            obj.createDesignVariable();
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
            x = ones(N+1,1);   
            s.type = 'AreaColumn';
            des = DesignVariable.create(s);
            des.update(x);
            obj.designVariable = des;  
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
            s.nElem          = obj.nElem;
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
            solution = IterativeProcessComputer(s);
            solution.compute();
            obj.designVariable = solution.designVariable;
        end     
        
    end
end