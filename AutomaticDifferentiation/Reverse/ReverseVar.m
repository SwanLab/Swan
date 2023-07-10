classdef ReverseVar < handle
    % It is just a test

    properties (Access = public)
        value
        children
        children2ndOrder
        grad_value
    end

    methods (Access = public)

        function obj = ReverseVar(val)
            obj.value = val;
        end

        function grad = computeGradient(obj)
            if isempty(obj.grad_value)
                obj.grad_value = 0;
                for i = 1:length(obj.children)
                    % disp('times & plus')
                    c = obj.children{i};
                    obj.grad_value = obj.grad_value + sum(c{1}.*c{2}.computeGradient());
                end
            end
            grad = obj.grad_value;
        end

        function z = mtimes(a,b)
            % disp('times')
            if ~isa(a,'ReverseVar') && ~isa(b,'ReverseVar')
                z = a*b;
            elseif ~isa(a,'ReverseVar')
                z = ReverseVar(b.value * a);
                b.children{end+1} = {a,z};
            elseif ~isa(b,'ReverseVar')
                z = ReverseVar(a.value * b);
                a.children{end+1} = {b,z};
            else
                z = ReverseVar(a.value * b.value);
                a.children{end+1} = {b.value, z, [0,1]};
                b.children{end+1} = {a.value, z, [0,1]};
            end
        end

        function z = plus(a,b)
            % disp('plus')
            if ~isa(a,'ReverseVar') && ~isa(b,'ReverseVar')
                z = a + b;
            elseif ~isa(a,'ReverseVar')
                z = ReverseVar(b.value + a);
                b.children{end+1} = {1,z};
            elseif ~isa(b,'ReverseVar')
                z = ReverseVar(a.value + b);
                a.children{end+1} = {1,z};
            else
                z = ReverseVar(a.value + b.value);
                a.children{end+1} = {1, z};
                b.children{end+1} = {1, z};
            end
        end

        function z = minus(a,b)
            % disp('minus')
            if ~isa(a,'ReverseVar') && ~isa(b,'ReverseVar')
                z = a - b;
            elseif ~isa(a,'ReverseVar')
                z = ReverseVar(b.value - a);
                b.children{end+1} = {-1,z};
            elseif ~isa(b,'ReverseVar')
                z = ReverseVar(a.value - b);
                a.children{end+1} = {1,z};
            else
                z = ReverseVar(a.value - b.value);
                a.children{end+1} = {1, z};
                b.children{end+1} = {-1, z};
            end
        end

        function z = mpower(a,b)
            % disp('power')
            if ~isa(a,'ReverseVar') && ~isa(b,'ReverseVar')
                z = a^b;
            elseif ~isa(a,'ReverseVar')
                z = ReverseVar(a^b.value);
                b.children{end+1} = {a^b.value*b.value*log(a),z};
            elseif ~isa(b,'ReverseVar')
                z = ReverseVar(a.value^b);
                a.children{end+1} = {b*a.value^(b-1),z};
            else
                z = ReverseVar(a.value^b.value);
                a.children{end+1} = {b.value*a.value^(b.value - 1), z};
                b.children{end+1} = {a.value^b.value*log(a.value), z};
            end
        end

        function z = sin(a)
            % disp('sin')
            if isa(a,'ReverseVar')
                z = ReverseVar(sin(a.value));
                a.children{end+1} = {cos(a.value),z};
            else
                z = sin(a);
            end
        end

        function z = cos(a)
            % disp('cos')
            if isa(a,'ReverseVar')
                z = ReverseVar(cos(a.value));
                a.children{end+1} = {-sin(a.value),z};
            else
                z = cos(a);
            end
        end

        function z = log(a)
            if isa(a,'ReverseVar')
                z = ReverseVar(log(a.value));
                a.children{end+1} = {1./a.value,z};
            else
                z = log(a);
            end
        end

    end

    methods (Access = public)



    end

end