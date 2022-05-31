classdef ShFunc_firstEigTopConst < ShapeFunctional
    
    properties (Access = private)
        msh
        dm
        eigModes
    end

    methods (Access = public)

        function obj = ShFunc_firstEigTopConst(cParams)
            obj.init(cParams)
            obj.msh = cParams.mesh;
            obj.createDim()
            obj.createEigModes();
        end

        function computeFunctionAndGradient(obj)
            obj.computeFunction();
            obj.computeGradient();
        end

        function t = getTitlesToPlot(obj)
            t{1} = 'First Eigen Value';
        end  

        function v = getVariablesToPlot(obj)
            v{1} = obj.value; % *obj.value0;
        end 
        
    end
% 
%     methods (Access = protected)
% 
%         function init(obj,cParams)
%             obj.createFilter(cParams);
%             % obj.createMsmoothAndDvolu(cParams);
%             obj.homogenizedVariablesComputer = cParams.homogVarComputer;
%             obj.designVariable = cParams.designVariable;
%             obj.target_parameters = cParams.targetParameters;
%             obj.nVariables = obj.designVariable.nVariables;     
%         end
% 
%     end

    methods (Access = private)

        function createDim(obj)
            s.type = 'Vector';
            s.ndimf = 2;
            s.name = 'x';
            s.fieldName = 'Disp';
            s.mesh = obj.msh;
            dims   = DimensionVariables.create(s);
            dims.compute();
            obj.dm = dims;
        end

        function createEigModes(obj)
            s.mesh           = obj.msh;
            s.designVariable = obj.designVariable;
            s.dim            = obj.dm;
            obj.eigModes = EigModesTopology(s);
        end

    end

    methods (Access = public)

        function computeFunction(obj)
            obj.value = obj.eigModes.provideFunction();
            % obj.normalizeFunction();
        end

        function computeGradient(obj)
            dfdx = obj.eigModes.provideDerivative();
            obj.gradient = dfdx';
            % obj.normalizeGradient();
        end
    end
    
end