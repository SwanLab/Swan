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

    % use dimension class for dim
    % length with geometry    
    % stiffnes and bending with LHSintegrator..
    % derivative "clean"/ "understand"    
    % MMa from Swan
    % Plot column area



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

            s.mesh = obj.mesh;


           quad = Quadrature.set(obj.mesh.type);
           quad.computeQuadrature('LINEAR');

            q   = quad;
            int = Interpolation.create(obj.mesh,'LINEAR');
            int.computeShapeDeriv(q.posgp);

            g = Geometry.create(s);   
            g.computeGeometry(q,int)    
            l = sum(g.dvolu,2);
        end

        function createDimensions(obj)
            s.mesh = obj.mesh;
            s.pdim = '1D'; 
            s.ngaus = 2;
            d = DimensionVariables(s);
            d.compute();
            d.ndof = 2*d.ndof; %%%%% Hay que arreglarlo !!!!!
            obj.dim = d;
        end

        function coord = createCoordinates(obj)
            nnode = obj.nElem + 1;
            x = [0;rand(nnode-2,1);1]*obj.columnLength;
            x = sort(x);
            % coord = x;
            coord = linspace(0,obj.columnLength,nnode)';
        end

        function Tnod = createConnectivity(obj)
            Tnod = zeros(obj.nElem,obj.nNodE);
            e = 1;
            for iElem = 1: obj.nElem
                Tnod(iElem,1) = e;
                e = e + 1;
                Tnod(iElem,2) = e;
            end
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
            obj.designVariable = solution.designVariable;
        end     
        
    end
end