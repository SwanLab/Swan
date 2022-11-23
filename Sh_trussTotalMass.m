classdef Sh_trussTotalMass < ShapeFunctional
    
    properties (Access = public)
        sectionAreaInfo
        barsLength
    end
    
    methods (Access = public)
        
        function obj = Sh_trussTotalMass(cParams)
            obj.init(cParams);
            % Load barres
            obj.homogenizedVariablesComputer.interpolation.loadFiles();
        end
        
        function computeFunctionAndGradient(obj)
            obj.nVariables = obj.designVariable.nVariables;
            obj.computeFunction();
            obj.computeGradient();
        end
        
        function computeFunction(obj)
            sect = obj.homogenizedVariablesComputer.interpolation.sectionArea;
            l    = obj.barsLength;
            obj.value = sect'*l;
        end
        
        function v = getVariablesToPlot(obj)
            v{1} = obj.value*obj.value0;
        end

        function t = getTitlesToPlot(obj)
            t{1} = 'Structure total mass';
        end
        
    end
    
    methods (Access = private)

        function computeGradient(obj)
            dSect = obj.homogenizedVariablesComputer.interpolation.lengthLessAreaGradient;
            l     = obj.barsLength;
            obj.gradient = dSect.*l;
        end

    end
    
end