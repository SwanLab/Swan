classdef ValDerForward

    properties

        val % function value

        der % derivative value or gradient vector

    end

    methods (Access = public)

        function obj = ValDerForward(a,b) %VALDER class constructor;

            obj.val = a; %given function value
            obj.der = b; %given derivative value or gradient vector

        end

        function vec = double(obj) %VALDER/DOUBLE Convert valder object to vector of doubles.

            vec = [ obj.val, obj.der ];

        end

        function h = plus(u,v)

            h = plusFun(u,v);

        end

        function h = minus(u,v)

            h = minusFun(u,v);

        end

        function h = uminus(u,v)

            h = uminusFun(u,v);

        end

        function h = mtimes(u,v)

            h = mtimesFun(u,v);

        end

        function h = mrdivide(u,v)

            h = mrdivideFun(u,v);

        end

        function h = mpower(u,v)

            h = mpowerFun(u,v);

        end

        function h = sin(u)

            h = sinFun(u);

        end

        function h = cos(u)

            h = cosFun(u);

        end

        function h = exp(u)

            h = expFun(u);

        end
    end
end