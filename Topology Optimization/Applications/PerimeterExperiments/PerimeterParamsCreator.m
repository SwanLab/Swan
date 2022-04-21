classdef PerimeterParamsCreator < handle
    
    properties (Access = public)
       perimeterParams 
    end
    
    properties (Access = private)
       inputFile
       mesh
       designVariable 
       epsilon
       scale
       type
       isRobinTermAdded
    end
    
    methods (Access = public)
        
        function obj = PerimeterParamsCreator(cParams)
            obj.init(cParams)
            obj.createPerimeterParams();
            obj.addFemParams();
            obj.addFilterParams();
            obj.addTargetParamters();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.inputFile        = cParams.inputFile;
            obj.mesh             = cParams.mesh;
            obj.designVariable   = cParams.designVariable;
            obj.epsilon          = cParams.epsilon;
            obj.scale            = cParams.scale;
            obj.type             = cParams.type;
            obj.isRobinTermAdded = cParams.isRobinTermAdded;
        end
        
        function createPerimeterParams(obj)
            obj.perimeterParams = SettingsShapeFunctional();
            obj.perimeterParams.type = 'perimeterInterior';
            obj.perimeterParams.designVariable = obj.designVariable;
            obj.perimeterParams.mesh = obj.mesh;
        end
           
        function addFemParams(obj)
           s = obj.createFemParams();
           obj.perimeterParams.femSettings = s;
        end
        
        function addFilterParams(obj)
            s = SettingsFilter();
            s.filterType =  'PDE';
            s.domainType =  'INTERIOR';
            s.femSettings = obj.createFemParams();
            s.designVarType = 'LevelSet';
            s.mesh = obj.mesh;
            obj.perimeterParams.filterParams = s;
        end
        
        function s = createFemParams(obj)
            s.fileName = obj.inputFile;
            s.scale    = obj.scale;
            s.mesh     = obj.mesh;
            s.isRobinTermAdded = obj.isRobinTermAdded;
        end
        
        function addTargetParamters(obj)
            s = TargetParameters();
            s.epsilon_perimeter = obj.epsilon;
            s.epsilon = obj.epsilon;
            obj.perimeterParams.targetParameters = s;
        end
        
    end
    
end