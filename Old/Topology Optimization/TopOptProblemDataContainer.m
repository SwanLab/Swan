classdef TopOptProblemDataContainer < AbstractSettings
    
    properties (Access = protected)
        defaultParamsName = ''
    end
    
    properties (Access = public)
        caseFileName
        femData
        
        costFunctions
        costWeights
        constraintFunctions
        nConstraints
    end
    
    methods (Access = public)
        
        function obj = TopOptProblemDataContainer(varargin)
            if nargin == 1
                obj.loadParams(varargin{1});
            end
            obj.init();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            s = obj.femData;
            obj.femData = FemDataContainer(s);
        end
        
        function readFemInputFile(obj)
            fileName = obj.femData.fileName;
            femReader = FemInputReader_GiD();
            s = femReader.read(fileName);
            
            obj.mesh  = s.mesh;
            obj.femData.scale = s.scale;
            obj.femData.dim  = s.pdim;
            obj.femData.type = s.ptype;
            obj.femData.nelem = s.mesh.nelem;
            obj.femData.bc.dirichlet = s.dirichlet;
            obj.femData.bc.pointload = s.pointload;
        end
        
    end
    
end