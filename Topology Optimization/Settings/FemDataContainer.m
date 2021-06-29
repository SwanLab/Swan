
classdef FemDataContainer < AbstractSettings
    
    properties (Access = protected)
        defaultParamsName = ''
    end
    
    properties (Access = public)
        fileName
        scale
        pdim
        ptype
        nelem
        bc
        mesh
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
            end
        end
        
        function readFemInputFile(obj)
            femReader = FemInputReader_GiD();
            s = femReader.read(obj.fileName);
            
            obj.mesh   = s.mesh;
            obj.scale  = s.scale;
            obj.pdim   = s.pdim;
            obj.ptype  = s.ptype;
            obj.nelem  = s.mesh.nelem;
            obj.bc.dirichlet = s.dirichlet;
            obj.bc.pointload = s.pointload;
        end
        
    end
    
end