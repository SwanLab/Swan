classdef ShFunc_firstEigTopConst < ShapeFunctional
    
    properties (Access = private)
        msh
        dm
        eigModes
        HomoVariablComputer
    end

    methods (Access = public)

        function obj = ShFunc_firstEigTopConst(cParams)
            obj.init(cParams)
            obj.msh = cParams.mesh;
            obj.createDim()
            obj.createEigModes();
            obj.HomoVarComputer();
        end

        function computeFunctionAndGradient(obj)
            obj.computeFunction();
            obj.computeGradient();
        end

        function t = getTitlesToPlot(obj)
            t{1} = 'First Eigenvalue Constraint';
        end  

        function v = getVariablesToPlot(obj)
            v{1} = obj.value; % *obj.value0;
        end 
        
    end

    methods (Access = protected)

        function HomoVarComputer(obj)
            obj.HomoVariablComputer = obj.homogenizedVariablesComputer;
        end

    end

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
            s.Msmooth        = obj.Msmooth;
            s.filter         = obj.filter;
            obj.eigModes = EigModesTopology(s);
        end

    end

    methods (Access = public)

        function computeFunction(obj)
            obj.va
            lue = obj.eigModes.provideFunction(obj.HomoVariablComputer);
        end

        function computeGradient(obj)
            dfdx = obj.eigModes.provideDerivative();
            obj.gradient = dfdx';
        end
    end
    
end