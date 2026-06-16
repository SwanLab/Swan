classdef MicroWithDensity < DesignVariable

    properties (Access = private)
        plotting
        plotter
        % funB
        % funRho
    end
    properties (Access = public)
        funB
        funRho
    end

    methods (Access = public)

        function obj = MicroWithDensity(cParams)
            obj.nVariables = 2;
            obj.initConcatenated(cParams);
            obj.funB   = cParams.funB;
            obj.funRho = cParams.fun;
            obj.createPlotter(cParams);
        end

        function fun = obtainDomainFunction(obj)
           
            fun{1} = obj.funB;
            fun{2} = obj.funRho;
        end

        function update(obj, value)
            
            n = size(obj.funB.fValues, 1);
            obj.funB.setFValues(value(1:n));
            obj.funRho.setFValues(value(n+1:end));
            obj.fun.setFValues(value);  
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

        function initConcatenated(obj, cParams)
            
            vals = [cParams.funB.fValues; cParams.fun.fValues];
            sFun        = cParams.fun.copy();
            sFun.setFValues(vals);
            cParams.fun = sFun;
            obj.init(cParams);   
        end

        function createPlotter(obj, cParams)
            obj.plotting = cParams.plotting;
            if obj.plotting
                obj.plotter = Plotter.create(obj);
            end
        end

    end

end