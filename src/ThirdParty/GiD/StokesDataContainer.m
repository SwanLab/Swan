classdef StokesDataContainer < AbstractSettings
    
    properties (Access = protected)
        defaultParamsName = ''
    end
    
    properties (Access = public)
        fileName
        scale
        state
        dtime
        finalTime
        dim
        type
        nelem
        bc
        mesh
        material
    end
    
    methods (Access = public)
        
        function obj = StokesDataContainer(varargin)
            if nargin == 1
                obj.loadParams(varargin{1}); %Fins aquí s'ha creat l'objecte que contè els paràmetres com la malla (però buits)
            end
            obj.init();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            if ~isempty(obj.fileName)
                obj.readFemInputFile();
                obj.createMaterial();
            end
        end
        
        function readFemInputFile(obj)
            femReader = FemInputReaderGiD();
            s = femReader.read(obj.fileName);
            
            obj.mesh      = s.mesh;
            obj.scale     = s.scale;
            obj.state     = s.state;
            obj.dtime     = s.dtime;
            obj.finalTime = s.ftime;
            obj.dim       = s.pdim;
            obj.type      = s.ptype;
            obj.nelem     = s.mesh.nelem;
            obj.bc.velocity = s.velocity;
            obj.bc.pressure = s.pressure;
            obj.bc.forcesFormula = s.forcesFormula;
            obj.bc.velocityBC    = s.velocityBC;
            obj.bc.dirichletFun = s.dirichletFun;
            obj.bc.pointloadFun = s.pointloadFun;
        end

        function createMaterial(obj)
            s.type = 'STOKES';
            s.nelem = obj.nelem;
            mat = Material.create(s);
            mat.compute();
            obj.material = mat;
        end

    end
    
end