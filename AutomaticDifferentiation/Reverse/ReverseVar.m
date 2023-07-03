classdef ReverseVar < handle
    % It is just a test

    properties (Access = public)
        value
        children
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

        function z = mtimes(obj,other)
            % disp('times')
            if isa(other,'ReverseVar')
                z = ReverseVar(obj.value * other.value);
                obj.children{end+1} = {other.value, z};
                other.children{end+1} = {obj.value, z};
            else
                z = ReverseVar(obj.value * other);
                obj.children{end+1} = {other,z};
            end
        end

        function z = plus(obj,other)
            % disp('plus')
            if isa(other,'ReverseVar')
                z = ReverseVar(obj.value + other.value);
                obj.children{end+1} = {1, z};
                other.children{end+1} = {1, z};
            else
                z = ReverseVar(obj.value + other);
                obj.children{end+1} = {1,z};
            end
        end

        function z = minus(obj,other)
            % disp('minus')
            if isa(other,'ReverseVar')
                z = ReverseVar(obj.value - other.value);
                obj.children{end+1} = {1, z};
                other.children{end+1} = {-1, z};
            else
                z = ReverseVar(obj.value - other);
                obj.children{end+1} = {1,z};
            end
        end

        function z = mpower(obj,other)
            % disp('power')
            if isa(other,'ReverseVar')
                z = ReverseVar(obj.value.^other.value);
                obj.children{end+1} = {other.value.*obj.value.^(other.value - 1), z};
                other.children{end+1} = {obj.value.^other.value.*log(obj.value), z};
            else
                z = ReverseVar(obj.value.^other);
                obj.children{end+1} = {other.*obj.value.^(other - 1),z};
            end
        end

        function z = sin(obj)
            % disp('sin')
            z = ReverseVar(sin(obj.value));
            obj.children{end+1} = {cos(obj.value),z};
        end

        function z = cos(obj)
            % disp('cos')
            z = ReverseVar(cos(obj.value));
            obj.children{end+1} = {-sin(obj.value),z};
        end

        function z = log(obj)
            z = ReverseVar(log(obj.value));
            obj.children{end+1} = {1./obj.value,z};
        end

    end

    methods (Access = public)

        

    end

end