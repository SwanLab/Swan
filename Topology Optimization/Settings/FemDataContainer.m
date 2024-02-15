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
            obj.init();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            if ~isempty(obj.fileName)
                obj.readFemInputFile();
                obj.getNgaus();
                obj.createMaterial();
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

        function createMaterial(obj)
            I = ones(obj.nelem,obj.ngaus);
            s.ptype = obj.type;
            s.pdim  = obj.dim;
            s.nelem = obj.nelem;
            s.mesh  = obj.mesh;
            s.kappa = .9107*I;
            s.mu    = .3446*I;
            mat = Material.create(s);
            mat.compute(s);
            obj.material = mat;
        end

        function getNgaus(obj)
            quad = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature('LINEAR');
            obj.ngaus = quad.ngaus;
        end

    end
    
end