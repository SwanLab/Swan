classdef FemDataContainer < AbstractSettings
    
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
        ngaus
        interpolationType
        solverType = 'REDUCED';
        solverMode = 'DISP';
        newBC
        boundaryConditions
    end
    
    methods (Access = public)
        
        function obj = FemDataContainer(varargin)
            if nargin == 1
                obj.loadParams(varargin{1});
            end
            obj.init(varargin{1});
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            if ~isempty(obj.fileName)
                obj.readFemInputFile();
                obj.getNgaus();
                obj.createMaterial(cParams);
                if strcmp(obj.scale, 'MICRO')
                    obj.solverMode = 'FLUC';
                end
            end
        end
        
        function readFemInputFile(obj)
            femReader = FemInputReader_GiD();
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

        function createMaterial(obj,cParams)
            E1        = 1;
            nu1       = 1/3;
            E         = AnalyticalFunction.create(@(x) E1*ones(size(squeeze(x(1,:,:)))),1,obj.mesh);
            nu        = AnalyticalFunction.create(@(x) nu1*ones(size(squeeze(x(1,:,:)))),1,obj.mesh);
            s.ptype   = obj.type;
            s.pdim    = obj.dim;
            s.nelem   = obj.nelem;
            s.mesh    = obj.mesh;
            s.young   = E;
            s.poisson = nu;
            if ~isfield(cParams,'type')
                cParams.type = 'ISOTROPIC';
            end
            s.type = cParams.type;
            s.ndim = obj.mesh.ndim;
            mat = Material.create(s);
            obj.material = mat;
        end

        function getNgaus(obj)
            quad = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature('LINEAR');
            obj.ngaus = quad.ngaus;
        end

    end
    
end