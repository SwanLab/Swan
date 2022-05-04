classdef EulerBeamOptimizer < handle

    properties (Access = public)
        columnMesh
    end
    
    properties (Access = protected)
        optimizerType
    end
    
    properties (Access = private)
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

%% TO DO LIST (29 april - 5 may)

    %%% Refactoring
    % delete elemental loop LHS Bending (DONE)
    % delete elemental loop LHS Stiffness (DONE)
    % delete LHS from TopOpt and transform the good one from Swan (DONE)

    %%% Next
    % 1. derivative "clean"/ "understand" (DONE)  
    % 2. MMa from Swan 
    % 3. postprocess with GiD; (with Alex)
    % 4. create shape functions cubic 1D (Interpolation) (finally???)
    
    %%% Future
    % Use shape functions for elemental stifness and bending
%%
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
            obj.createDimensions();
            obj.createDesignVariable();
            obj.createBoundaryConditions()
            obj.createMMA();
            obj.computeIterativeProcess()
            obj.PostProcess();
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
            s.coord  = obj.createCoordinates();
            s.connec = obj.createConnectivity();
            s.type = 'LINE';
            m = Mesh(s);
            obj.mesh = m;
        end

        function coord = createCoordinates(obj)
            nnod = obj.nElem + 1;
             x = [0;rand(nnod-2,1);1]*obj.columnLength;
             x = sort(x);
             % coord = x;
            coord = linspace(0,obj.columnLength,nnod)';
        end

        function Tnod = createConnectivity(obj)
            nNode = 2;
            Tnod = zeros(obj.nElem,nNode);
            e = 1;
            for iElem = 1: obj.nElem
                Tnod(iElem,1) = e;
                e = e + 1;
                Tnod(iElem,2) = e;
            end
        end

        function createDimensions(obj)
            s.mesh = obj.mesh;
            s.pdim = '2D';
            s.ngaus = 2;
            d = DimensionVariables(s);
            d.compute();
            obj.dim = d;
        end

        function createDesignVariable(obj)
            s.type  = 'AreaColumn';
            s.mesh  = obj.mesh;
            des = DesignVariable.create(s);
            obj.designVariable = des;  
        end

        function createBoundaryConditions(obj)
            d = obj.dim;
            fixnodes = union([1,2], [d.ndof-1,d.ndof]);
            nodes = 1:d.ndof;
            free  = setdiff(nodes,fixnodes);
            obj.freeNodes = free;
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
            s.mmaVarComputer = obj.mmaVarComputer;
            s.freeNodes      = obj.freeNodes;
            s.nConstraints   = obj.nConstraints;
            s.mesh           = obj.mesh;
            s.nValues        = obj.nValues;
            s.dim            = obj.dim;
            s.youngModulus   = obj.youngModulus;
            s.inertiaMoment  = obj.inertiaMoment;
            s.maxIter        = obj.maxIter;
            s.optimizerType  = obj.optimizerType;
            s.designVariable = obj.designVariable;
            solution = IterativeProcessComputer(s);
            solution.compute();
        end

        function PostProcess(obj)
            s.designVariable = obj.designVariable;
            s.dim            = obj.dim;
            s.mesh           = obj.mesh;
            s.scale          = 0.3;
            post = PostProcessColumn(s);
            post.plotColumn();
            obj.columnMesh = post.m;
        end
        
    end
end