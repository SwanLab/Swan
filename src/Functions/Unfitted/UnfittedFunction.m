classdef UnfittedFunction < BaseFunction

    properties (Access = public)
        unfittedMesh
        innerMeshFunction
        innerCutMeshFunction
    end

    properties (Access = private)
        fun
    end

    methods (Static, Access = public)

        function obj = create(uMesh,fun)
            s.uMesh = uMesh;
            s.fun   = fun;
            obj     = UnfittedFunction(s);
        end

    end

    methods (Access = public)

        function obj = UnfittedFunction(cParams)
            obj.init(cParams);
            obj.computeUnfittedMeshFunction();
            obj.mesh = obj.unfittedMesh.backgroundMesh;
        end

        function res = times(obj1,obj2)
            res = copy(obj1);
            switch class(obj2)
                case 'UnfittedFunction'
                    res.fun = res.fun.*obj2.fun;
                    res.innerMeshFunction = res.innerMeshFunction.*obj2.innerMeshFunction;
                    res.innerCutMeshFunction = res.innerCutMeshFunction.*obj2.innerCutMeshFunction;
                    res.updateNDimF(obj2);
                case {'LagrangianFunction','AnalyticalFunction'}
                    res.fun = res.fun.*obj2;
                    f2      = obj1.createNew(obj2);
                    res.innerMeshFunction = res.innerMeshFunction.*f2.innerMeshFunction;
                    res.innerCutMeshFunction = res.innerCutMeshFunction.*f2.innerCutMeshFunction;
                    res.updateNDimF(f2);
                case 'double'
                    res.fun = res.fun.*obj2;
                    res.innerMeshFunction = res.innerMeshFunction.*obj2;
                    res.innerCutMeshFunction = res.innerCutMeshFunction.*obj2;
                case 'DomainFunction'
                    error('Can not be a domain function');
            end
        end

        function res = power(obj,p)
            if p == 0
                res = ConstantFunction.create(1,obj.fun.mesh);
            else
                res = copy(obj);
                res.fun = res.fun.^p;
                res.innerMeshFunction = res.innerMeshFunction.^p;
                res.innerCutMeshFunction = res.innerCutMeshFunction.^p;
            end
        end

        function res = plus(obj1,obj2)
            res = copy(obj1);
            switch class(obj2)
                case 'UnfittedFunction'
                    res.fun = res.fun + obj2.fun;
                    res.innerMeshFunction = res.innerMeshFunction + obj2.innerMeshFunction;
                    res.innerCutMeshFunction = res.innerCutMeshFunction + obj2.innerCutMeshFunction;
            end
        end

        function res = minus(obj1,obj2)
            res = copy(obj1);
            switch class(obj2)
                case 'UnfittedFunction'
                    res.fun = res.fun - obj2.fun;
                    res.innerMeshFunction = res.innerMeshFunction - obj2.innerMeshFunction;
                    res.innerCutMeshFunction = res.innerCutMeshFunction - obj2.innerCutMeshFunction;
            end
        end

        function res = rdivide(obj1,obj2)
            res = copy(obj1);
            switch class(obj2)
                case 'UnfittedFunction'
                    res.fun = res.fun./obj2.fun;
                    res.innerMeshFunction = res.innerMeshFunction./obj2.innerMeshFunction;
                    res.innerCutMeshFunction = res.innerCutMeshFunction./obj2.innerCutMeshFunction;
                case {'LagrangianFunction','AnalyticalFunction'}
                    res.fun = res.fun./obj2;
                    f2      = obj1.createNew(obj2);
                    res.innerMeshFunction = res.innerMeshFunction./f2.innerMeshFunction;
                    res.innerCutMeshFunction = res.innerCutMeshFunction./f2.innerCutMeshFunction;
                case 'double'
                    res.fun = res.fun./obj2;
                    res.innerMeshFunction = res.innerMeshFunction./obj2;
                    res.innerCutMeshFunction = res.innerCutMeshFunction./obj2;
                case 'DomainFunction'
                    error('Can not be a domain function');
            end
        end

    end

    methods (Access = private)
        function init(obj,cParams)
            obj.unfittedMesh = cParams.uMesh;
            obj.fun          = cParams.fun;
            obj.ndimf        = cParams.fun.ndimf;
        end

        function computeUnfittedMeshFunction(obj)
            uMeshFun = obj.unfittedMesh.obtainFunctionAtUnfittedMesh(obj.fun);
            obj.innerMeshFunction    = uMeshFun.innerMeshFunction;
            obj.innerCutMeshFunction = uMeshFun.innerCutMeshFunction;
        end

        function updateNDimF(obj,f2)
            obj.ndimf = max(obj.ndimf,f2.ndimf);
        end

        function f = createNew(obj,fun)
            s.uMesh = obj.unfittedMesh;
            s.fun   = fun;
            f       = UnfittedFunction(s);
        end
    end

    methods (Access = protected)

        function evaluateNew(obj,xV)

        end

    end
end