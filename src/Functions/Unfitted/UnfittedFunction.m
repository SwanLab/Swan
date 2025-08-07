classdef UnfittedFunction < BaseFunction

    properties (Access = public)
        unfittedMesh
        innerMeshFunction
        innerCutMeshFunction
    end

    properties (Access = private)
        fun
    end

    methods (Access = public)

        function obj = UnfittedFunction(cParams)
            obj.init(cParams);
            obj.computeUnfittedMeshFunction();
            obj.mesh = obj.unfittedMesh.backgroundMesh;
        end

        function res = times(obj1,obj2)
            res     = copy(obj1);
            switch class(obj2)
                case 'UnfittedFunction'
                    res.fun = res.fun.*obj2.fun;
                otherwise
                    res.fun = res.fun.*obj2;
            end
            res.computeUnfittedMeshFunction();
        end

        function res = power(obj,p)
            if p == 0
                res = ConstantFunction.create(1,obj.fun.mesh);
            else
                res = copy(obj);
                res.fun = res.fun.^p;
                res.computeUnfittedMeshFunction();
            end
        end

        function res = plus(obj1,obj2)
            res     = copy(obj1);
            switch class(obj2)
                case 'UnfittedFunction'
                    res.fun = res.fun + obj2.fun;
                    res.computeUnfittedMeshFunction();
            end
        end

        function res = minus(obj1,obj2)
            res     = copy(obj1);
            switch class(obj2)
                case 'UnfittedFunction'
                    res.fun = res.fun - obj2.fun;
                    res.computeUnfittedMeshFunction();
            end
        end

        function res = rdivide(obj1,obj2)
            res     = copy(obj1);
            switch class(obj2)
                case 'UnfittedFunction'
                    res = res.fun./obj2.fun;
                otherwise
                    res.fun = res.fun./obj2;
                    res.computeUnfittedMeshFunction();
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
    end

    methods (Access = protected)

        function evaluateNew(obj,xV)

        end

    end
end