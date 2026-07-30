classdef FemDataContainer < handle
    
    properties (Access = protected)
        defaultParamsName = ''
    end
    
    properties (Access = public)
        fileName
        scale
        dim
        type
        nelem
        bc
        mesh
        material
        interpolationType
        solverType = 'REDUCED';
        solverMode = 'DISP';
        solverCase = DirectSolver();
        newBC
        boundaryConditions
    end
    
    methods (Access = public)
        
        function obj = FemDataContainer(varargin)
            obj.fileName = varargin{1}.fileName;
            obj.init();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            if ~isempty(obj.fileName)
                obj.readFemInputFile();
                obj.createMaterial();
                if strcmp(obj.scale, 'MICRO')
                    obj.solverMode = 'FLUC';
                end
            end
        end
        
        function readFemInputFile(obj)
            femReader = FemInputReaderGiD();
            s = femReader.read(obj.fileName);

            obj.mesh   = s.mesh;
            obj.scale  = s.scale;
            obj.dim   = s.pdim;
            obj.type  = s.ptype;
            obj.nelem  = s.mesh.nelem;
            obj.bc.dirichlet = s.dirichlet;
            obj.bc.pointload = s.pointload;
            obj.interpolationType = 'LINEAR';
            obj.newBC.dirichletFun = s.dirichletFun;
            obj.newBC.pointloadFun = s.pointloadFun;
            obj.newBC.periodicFun  = s.periodicFun;
            obj.boundaryConditions = BoundaryConditions(s);
        end

        function createMaterial(obj)
            N      = obj.mesh.ndim;
            E      = ConstantFunction.create(1,obj.mesh);
            nu     = ConstantFunction.create(1/3,obj.mesh);
            mu     = E./(2*(1+nu));  mu = Expand(mu,4);
            lambda = LameLambda(E,nu,N);  lambda = Expand(lambda,4);
            I      = ConstantFunction.create(eye4D(N),obj.mesh);
            IxI    = ConstantFunction.create(kronEye(N),obj.mesh);
            obj.material = 2*mu.*I + lambda.*IxI;
        end

    end
    
end