classdef MaterialDesignExperiment < handle
    
   properties (Access = public)
        
    end
    
    properties (Access = private)
        
    end
    
    properties (Access = private)
        fileNames
        topOptSet
        topOptProblem
    end
    
    methods (Access = public)
        
        function obj = MaterialDesignExperiment()
            obj.init();
            for icases = 1:numel(obj.fileNames)
                obj.createSettings(icases);
                obj.solveProblem();
                close all
            end
            
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.fileNames = {...
                             'HorizontalMaterialDesign';
                          %   'CompositeMaterialDesignTriDensityP1';
                           %  'CompositeMaterialDesignTriDensityPDE';
                           %  'CompositeMaterialDesignQuadDensityP1';
                            %  'CompositeMaterialDesignQuadDensityPDE';
                            % 'CompositeMaterialDesignTriLevelSetP1';
                            % 'CompositeMaterialDesignTriLevelSetPDE';
                           % 'CompositeMaterialDesignQuadLevelSetP1';
                           %  'CompositeMaterialDesignQuadLevelSetPDE';
                             };
        end
        
        
        function createSettings(obj,icases)
            s = SettingsTopOptProblem(obj.fileNames{icases});
            obj.topOptSet = s;
        end
        
        function solveProblem(obj)
            obj.topOptProblem = TopOpt_Problem(obj.topOptSet);
            obj.topOptProblem.computeVariables();
            obj.topOptProblem.postProcess();
        end
    end
    
end