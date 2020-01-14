classdef PerimeterParamsCreator < handle
    
    properties (Access = public)
       filterParams 
    end
    
    properties (Access = private)
       inputFile
       mesh
       designVariable 
       epsilon
    end
    
    methods (Access = public)
        
        function obj = PerimeterParamsCreator(cParams)
            obj.init(cParams)
            obj.createFilterParams()
            obj.createFemParams();
            obj.createPerimeterParams();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.inputFile      = cParams.inputFile;
            obj.mesh           = cParams.mesh;
            obj.designVariable = cParams.designVariable;
            obj.epsilon        = cParams.epsilon;
        end
        
        function createPerimeterParams(obj)
            obj.perimeterParams = SettingsShapeFunctional();      
            obj.perimeterParams.type = 'perimeterInterior';            
            obj.perimeterParams.designVariable   = obj.designVariable;            
        end
           
        function createFilterParams(obj)
            s = SettingsFilter();
            s.filterType =  'PDE';
            s.domainType =  'INTERIOR';
            s.designVar  =  obj.designVariable;
            s.quadratureOrder =  'LINEAR';
            s.femSettings = obj.createFemParams();
            obj.filterParams.filterParams = s;            
        end        
        
        function s = createFemParams(obj)
            s.fileName = obj.inputFile;
            s.scale    = 'MACRO';
            s.mesh     = obj.mesh;
            s.isRobinTermAdded = false;
        end        
        
        function createTargetParamters(obj)
            s = TargetParameters();
            s.epsilon_perimeter = obj.epsilon;
            s.epsilon = obj.epsilon;
            obj.filterParams.targetParameters = s;            
        end
        
        function createHomogenizedVarComputer(obj)
            s = SettingsHomogenizedVarComputer;
            obj.filterParams.homogVarComputer = s;
        end        
                        
    end
    
    
    
    
    
    
    
    
end