classdef ValGradReverse

    properties

        val

        Grad

        adj

    end

    methods

        function obj = ValGradReverse(value)

            obj.val = value;
            obj.Grad = 0;
            obj.adj = [1];

        end

        function h = plus(x, y)

            h = ValGradReverse(x.val + y.val);
            x.adj = [x.adj; {1, h}];
            y.adj = [y.adj; {1, h}];

        end


        function h = mtimes(x, y)

            h = ValGradReverse(x.val * y.val);
            x.adj = [x.adj; {y.val, h}];
            y.adj = [y.adj; {x.val, h}];

        end

        function h = mpower(u,v)

            %VALGrad/MPOWER overloads ^ with at least one valGrad
            if ~isa(u,'ValGradReverse') %u is a scalar
                h = ValGradReverse(u^v.val); %, u^v.val*log(u)*v.Grad);
            elseif ~isa(v,'ValGradReverse') %v is a scalar
                h = ValGradReverse(u.val^v); %, v*u.val^(v-1)*u.Grad);
            else
                h = exp(v*log(u)); %call overloaded log, * and exp
            end

        end


        function h = cos(x)

            h = ValGradReverse(cos(x.val));

        end


        function backward(obj, gradient)

            obj.Grad = obj.Grad + gradient;

            for i = 1:length(obj.adj)

                weight = obj.adj{i}{1};
                parent = obj.adj{i}{2};
                parent.backward(gradient * weight);

            end
        end
    end
end