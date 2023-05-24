classdef ValGradForward

    properties

        val % function value

        grad % Gradient value or gradient vector

    end

    methods (Access = public)

        function obj = ValGradForward(a,b) %ValGrad class constructor;
            if nargin == 0
            else
                if isequal(size(a),[3,3])
                     error
                end
                obj.val = a; %given function value
                obj.grad = b; %given Gradient value or gradient vector
            end

        end

        function vec = double(obj) %Convert ValGrad object to vector of doubles.

            if isequal(size(obj.val),[3,3]) || isequal(size(obj.grad),[3,3])
                vec = obj.val;
                vec(1:3,4:6) = obj.grad;

            else

                vec = [ obj.val, obj.grad ];
            end

        end

        function h = plus(u,v)

            h = plusFun(u,v);

        end

        function h = minus(u,v)

            h = minusFun(u,v);

        end

        function h = uminus(u)

            h = uminusFun(u);

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

        function h = log(u)

            h = logFun(u);

        end

        function h = sqrt(u)

            h = sqrtFun(u);

        end
    end
end