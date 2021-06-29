classdef IncrementalScheme < handle
    
    properties (GetAccess = public, SetAccess = private)
        iStep
        nSteps
        targetParams
    end
    
    properties (Access = private)
        targetParamsManager
        shallDisplayStep
    end
    
    methods (Access = public)
        
        function obj = IncrementalScheme(cParams)
            obj.init(cParams);
            obj.createTargetParams(cParams);
        end
        
        function next(obj)
            obj.incrementStep();
            obj.updateTargetParams();
        end
        
        function display(obj)
            disp(['Incremental Scheme - Step: ',int2str(obj.iStep),' of ',int2str(obj.nSteps)]);
        end
        
        function itDoes = hasNext(obj)
            if obj.iStep < obj.nSteps
                itDoes = true;
            else
                itDoes = false;
            end
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.iStep = 0;
            obj.nSteps = cParams.nSteps;
            obj.shallDisplayStep = cParams.shallPrintIncremental;
        end
        
        function createTargetParams(obj,cParams)
            s = cParams.targetParamsSettings;
            obj.targetParamsManager = TargetParamsManager(s);
            obj.targetParams = obj.targetParamsManager.targetParams;
        end
        
        function updateTargetParams(obj)
            obj.targetParamsManager.update(obj.iStep);
        end
        
        function incrementStep(obj)
            obj.iStep = obj.iStep + 1;
            if obj.shallDisplayStep
                obj.display();
            end
        end
        
    end
    
end