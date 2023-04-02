classdef ValDerForward

    properties

        val % function value

        der % derivative value or gradient vector

    end

    methods

        function obj = ValDerForward(a,b)

            %VALDER class constructor;

            obj.val = a; %given function value
            obj.der = b; %given derivative value or gradient vector

        end

        function vec = double(obj)

            %VALDER/DOUBLE Convert valder object to vector of doubles.

            vec = [ obj.val, obj.der ];

        end

        function h = plus(u,v)

             if ~isa(u,'ValDerForward') %u is a scalar
                h = ValDerForward(u+v.val, v.der);
            elseif ~isa(v,'ValDerForward') %v is a scalar
                h = ValDerForward(v+u.val, u.der);
            else
                h = ValDerForward(u.val+v.val, u.der+v.der);
             end

        end

        function h = minus(u,v)

            h = ValDerForward(u.val - v.val, u.der - v.der);

        end

        function h = mtimes(u,v)

            %VALDER/MTIMES overloads * with at least one valder

            if ~isa(u,'ValDerForward') %u is a scalar
                h = ValDerForward(u*v.val, u*v.der);
            elseif ~isa(v,'ValDerForward') %v is a scalar
                h = ValDerForward(v*u.val, v*u.der);
            else
                h = ValDerForward(u.val*v.val, u.der*v.val + u.val*v.der);
            end

        end

        function h = mpower(u,v)

            %VALDER/MPOWER overloads ^ with at least one valder
            if ~isa(u,'ValDerForward') %u is a scalar
                h = ValDerForward(u^v.val, u^v.val*log(u)*v.der);
            elseif ~isa(v,'ValDerForward') %v is a scalar
                h = ValDerForward(u.val^v, v*u.val^(v-1)*u.der);
            else
                h = exp(v*log(u)); %call overloaded log, * and exp
            end

        end

        function h = sin(u)

            h = ValDerForward(sin(u.val), cos(u.val)*u.der);

        end

        function h = cos(u)

            h = ValDerForward(cos(u.val), -sin(u.val)*u.der);

        end
    end
end