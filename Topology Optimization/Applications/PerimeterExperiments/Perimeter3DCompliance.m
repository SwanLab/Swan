classdef Perimeter3DCompliance < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        fileName
        topOptSet
        topOptProblem
    end
    
    methods (Access = public)
        
        function obj = Perimeter3DCompliance()
            obj.init();              
            obj.createTopOptSettings();            
            obj.solveProblem();
        end        
        
    end
    
    methods (Access = private)
        
        function init(obj)           
          %  obj.fileName = 'CantileverTetra';
            obj.fileName = 'CantileverFlat';
        end
                
        function createTopOptSettings(obj)
          t = SettingsTopOptProblem(obj.fileName);            
          obj.topOptSet = t;              
        end
        
        function solveProblem(obj)
            obj.topOptProblem = TopOpt_Problem(obj.topOptSet);
            obj.topOptProblem.computeVariables();
            obj.topOptProblem.postProcess();
        end        
        
    end
    
end