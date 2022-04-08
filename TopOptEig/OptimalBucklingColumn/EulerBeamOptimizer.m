classdef EulerBeamOptimizer < handle
    
    properties (Access = protected)
        optimizerType
    end
    
    properties (Access = private)
        nDim
        nElem
        nConstraints
        columnLength
        nNodes
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
       nDofN
       Tnod
       nNodE
   end

    % Mesh (lenght only to create mesh, delete elsewhere)  (DONE)
    % Kelem  + assembly    (DONE)
    % plot modes getting displacement (DONE)
    % Solve for a non-structured mesh (DONE)


    properties (Access = private)
        designVariable
        freeNodes
        nIter
        mmaVarComputer
    end
     
    methods (Access = public)
        
        function obj = EulerBeamOptimizer()
            obj.init()
            obj.createMesh();
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
            obj.columnLength  = 1; 
            obj.nValues       = obj.nElem+1;
            obj.youngModulus  = 1;
            obj.inertiaMoment = 1;  
            obj.minThick      = 0.25;
            obj.maxThick      = 10;
            obj.optimizerType = 'MMA';
            obj.maxIter       = 1000;
        end

        function createMesh(obj)
            obj.dim.nNod  = obj.nElem + 1;                  % num of nodes
            obj.mesh      = linspace(0,obj.columnLength,obj.dim.nNod)';
            obj.dim.nDim  = size(obj.mesh,2);               % dimension of the problem
            obj.dim.nNodE = 2;                              % num of node per element
            obj.Tnod      = obj.nodesConnectivityMatrix();
            obj.dim.nDofN = 2;                              % num of DOFs per element
            obj.dim.nDof  = obj.dim.nNod*obj.dim.nDofN;     % num of DOFs
            
        end

        function Tnod = nodesConnectivityMatrix(obj) % better to be an input
            Tnod = zeros(obj.nElem,obj.nNodE);
            e = 1;
            for iElem = 1: obj.nElem
                Tnod(iElem,1) = e;
                e = e + 1;
                Tnod(iElem,2) = e;
            end
        end

        function createDesignVariable(obj)
            N = obj.nElem;
            s.type  = 'AreaColumn';
            s.nElem = obj.nElem;
            s.mesh  = obj.mesh;
            des = DesignVariable.create(s);
            x0 = ones(N+1,1);               
            des.update(x0);
            obj.designVariable = des;  
        end

        function createBoundaryConditions(obj)
            d = obj.dim;
            fixnodes = union([1,2], [d.nDof-1,d.nDof]);
            nodes = 1:d.nDof;
            free  = setdiff(nodes,fixnodes);
            obj.freeNodes = free;
        end        

        function obj = createMMA(obj)
            s.nElem         = obj.nElem;
            s.nConstraints  = obj.nConstraints;
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
            s.mesh           = obj.mesh;
            s.Tnod           = obj.Tnod;
            s.nValues        = obj.nValues;
            s.dim            = obj.dim;
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