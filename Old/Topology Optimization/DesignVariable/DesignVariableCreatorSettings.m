classdef DesignVariableCreatorSettings < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        
    end
    
    properties (Access = private)
       mesh
       inputFile
       scale
       type
    end
    
    methods (Access = public)
        
        function obj = DesignVariableCreatorSettings(cParams)
            obj.init(cParams)
        end
        
        function s = create(obj)
            s = SettingsDesignVariable();
            s.mesh                    = Mesh_Total(obj.mesh);
            s.type                    = obj.type;
            s.initialCase             = 'full';
            s.creatorSettings         = obj.createLevelSetParams();
            s.scalarProductSettings   = obj.createScalarProductParams();
            s.femData                 = obj.createFemContainerData();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh      = cParams.mesh;
            obj.inputFile = cParams.inputFile;
            obj.scale     = cParams.scale;
            obj.type      = cParams.type;
        end
        
        function s = createFemContainerData(obj)
            s = FemDataContainer;
            s.fileName = obj.inputFile;
            s.scale    = obj.scale;
            s.pdim     = '2D';
            s.ptype    = 'ELASTIC';
            s.nelem    = size(obj.mesh.connec,1);
            s.bc       = [];
        end
        
        function s = createScalarProductParams(obj)
            s.scalarProductSettings.femSettings = [];
            s.epsilon = [];
            s.mesh = obj.mesh;
        end
        
        function s = createLevelSetParams(obj)
            ss.type = 'full';
            s = SettingsLevelSetCreator;
            s = s.create(ss);
        end
    end
    
end