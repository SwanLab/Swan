classdef Sh_trussTotalMass < ShapeFunctional
    
    properties (Access = public)
        sectionAreaInfo
        barsLength
    end
    
    methods (Access = public)
        
        function obj = Sh_trussTotalMass(cParams)
            cParams.filterParams.quadratureOrder = 'LINEAR';
            obj.init(cParams);
%             obj.geometricVolume = sum(obj.dvolu(:));
        end
        
        function computeFunctionAndGradient(obj)
            obj.nVariables = obj.designVariable.nVariables;
            obj.updateHomogenizedMaterialProperties();
            obj.computeFunction();
            obj.computeGradient();
        end
        
        function computeFunction(obj)
            sect = obj.homogenizedVariablesComputer.rho; % Això ni flis
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
            dSect = obj.homogenizedVariablesComputer.drho; % Això ni flis
            l     = obj.barsLength;
            obj.gradient = dSect.*l;
        end
        
%         function updateHomogenizedMaterialProperties(obj)
%             nx = length(obj.designVariable.value)/obj.designVariable.nVariables;
%             x  = obj.designVariable.value;
%             xf = cell(obj.designVariable.nVariables,1);
%             for ivar = 1:obj.nVariables
%                 i0 = nx*(ivar-1) + 1;
%                 iF = nx*ivar;
%                 xs = x(i0:iF);
%                 xf{ivar} = obj.filter.getP0fromP1(xs);
%             end
%             obj.homogenizedVariablesComputer.computeDensity(xf);
%         end

    end
    
end