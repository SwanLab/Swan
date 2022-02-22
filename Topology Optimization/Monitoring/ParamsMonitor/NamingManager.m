classdef NamingManager < handle
    
    properties (Access = private)
        costFuncNames
        costWeights
        constraintNames
        convVarsNames
    end
    
    methods (Access = public)
        
        function obj = NamingManager(cParams)
            obj.costFuncNames   = cParams.costFuncNames;
            obj.costWeights     = cParams.costWeights;
            obj.constraintNames = cParams.constraintFuncs;
            obj.convVarsNames   = ConvergenceVarsDispatcher.dispatchNames(cParams.optimizerNames);
        end
        
        function name = getCostFuncFigureTitle(obj,i)
            name = obj.costFuncNames{i};
            weight = obj.costWeights(i);
            
            if isempty(weight)
                weightDesc = '(wt. 1.0)';
            else
                weightDesc = sprintf('(wt. %.2f)',weight);
            end
            
            name = [obj.setCase(name) ' ' weightDesc];
        end
        
        function name = getConstraintFigureTitle(obj,i)
            name = obj.constraintNames{i};
            name = ['Cstr ' num2str(i) ': ' obj.setCase(name)];
        end
        
        function name = getLambdaFigureTitle(obj,i)
            name = obj.constraintNames{i};
            name = ['\lambda' obj.subindex(obj.setCase(name))];
        end
        
        function name = getConvVarFigureTitle(obj,i)
            name = obj.convVarsNames{i};
            name = ['Conv Criteria ' num2str(i) ': ' name];
        end
        
        function n = getNumberOfConvergenceVariables(obj)
            n = numel(obj.convVarsNames);
        end
        
        function name = getConvVarNames(obj,i)
            name = obj.convVarsNames{i};
        end
        
    end
    
    methods (Access = public, Static)
        
        function name = getFrameTitle(iStep,nStep,it)
            name = sprintf('Monitoring - Inc. Step: %.0f/%.0f Iteration: %.0f',iStep,nStep,it);
        end
        
    end
    
    methods (Access = private,Static)
        
        function str = setCase(str)
            str = [upper(str(1)) str(2:end)];
        end
        
        function str = subindex(str_)
            str = [];
            for i = 1:length(str_)
                str(end+1:end+2) = ['_' str_(i)];
            end
        end
    end
    
end

