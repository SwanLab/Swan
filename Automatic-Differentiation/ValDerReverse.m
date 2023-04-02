classdef ValDerReverse

    properties

        val

        der

        adj

    end

    methods

        function obj = ValDerReverse(value)

            obj.val = value;
            obj.der = 0;
            obj.adj = [1];

        end

        function h = plus(x, y)

            h = ValDerReverse(x.val + y.val);
            x.adj = [x.adj; {1, h}];
            y.adj = [y.adj; {1, h}];

        end


        function h = mtimes(x, y)

            h = ValDerReverse(x.val * y.val);
            x.adj = [x.adj; {y.val, h}];
            y.adj = [y.adj; {x.val, h}];

        end

        function h = mpower(u,v)

            %VALDER/MPOWER overloads ^ with at least one valder
            if ~isa(u,'ValDerReverse') %u is a scalar
                h = ValDerReverse(u^v.val); %, u^v.val*log(u)*v.der);
            elseif ~isa(v,'ValDerReverse') %v is a scalar
                h = ValDerReverse(u.val^v); %, v*u.val^(v-1)*u.der);
            else
                h = exp(v*log(u)); %call overloaded log, * and exp
            end

        end


        function h = cos(x)

            h = ValDerReverse(cos(x.val));

        end


        function backward(obj, gradient)

            obj.der = obj.der + gradient;

            for i = 1:length(obj.adj)

                weight = obj.adj{i}{1};
                parent = obj.adj{i}{2};
                parent.backward(gradient * weight);

            end
        end
    end
end