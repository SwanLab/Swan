classdef MicroWithDensity < DesignVariable

    properties (Access = private)
        plotting
        plotter
        funB
        funRho
        fValuesOld
        funConcatenated  % proxy com 2n valores para o ProjectedGradient
    end

    methods (Access = public)

        function obj = MicroWithDensity(cParams)
            obj.nVariables = 2;
            obj.init(cParams);
            obj.createPlotter(cParams);
            obj.funB   = cParams.funB;
            obj.funRho = cParams.fun;
            obj.createConcatenatedFun();
        end

        function fun = obtainDomainFunction(obj)
            fun{1} = obj.funB;
            fun{2} = obj.funRho;
        end

        function update(obj, value)
            n = size(obj.funB.fValues, 1);
            obj.funB.setFValues(value(1:n));
            obj.funRho.setFValues(value(n+1:end));
            obj.funConcatenated.setFValues(value);
            obj.fun = obj.funConcatenated;
        end

        function updateOld(obj)
            obj.fValuesOld = [obj.funB.fValues; obj.funRho.fValues];
        end

        function recoverOld(obj)
            n = length(obj.fValuesOld)/2;
            obj.funB.setFValues(obj.fValuesOld(1:n));
            obj.funRho.setFValues(obj.fValuesOld(n+1:end));
            obj.funConcatenated.setFValues(obj.fValuesOld);
            obj.fun = obj.funConcatenated;
        end

        function norm = computeL2normIncrement(obj)
            xNew = [obj.funB.fValues; obj.funRho.fValues];
            xOld = obj.fValuesOld;
            diff = xNew - xOld;
            norm = sqrt(sum(diff.^2)) / (sqrt(sum(xOld.^2)) + 1e-12);
        end

        function fixedDofs = getFixedDofs(obj)
            fixedDofs = obj.isFixed;
        end

        function plot(obj)
            if obj.plotting
                obj.plotter.plot();
            end
        end

    end

    methods (Access = private)

        function createConcatenatedFun(obj)
            % Cria um LagrangianFunction com 2n valores
            % copiando funB e concatenando com funRho
            obj.funConcatenated = obj.funB.copy();
            vals = [obj.funB.fValues; obj.funRho.fValues];
            obj.funConcatenated.setFValues(vals);
            obj.fun = obj.funConcatenated;
        end

        function createPlotter(obj, cParams)
            obj.plotting = cParams.plotting;
            if obj.plotting
                obj.plotter = Plotter.create(obj);
            end
        end

    end
end