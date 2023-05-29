classdef SymmetricLHSComputer < LHSComputer

    properties (Access = public)
        LHS
    end
    properties (Access = private)
    end

    methods (Access = public)
        function obj = SymmetricLHSComputer(cParams)
            obj.init(cParams);
        end
        function compute(obj)
            obj.computeLHS();
        end
    end

    methods (Access = private)
        function computeLHS(obj)
            sizeH = size(obj.H,1);
            sizeg = obj.m;
            lhs = zeros(sizeH+sizeg);
            lhs(1:sizeH,1:sizeH) = obj.H;
            lhs(1:sizeH,sizeH+1:sizeH+sizeg) = obj.constraint.gradient;
            lhs(sizeH+1:sizeH+sizeg,1:sizeH) = obj.constraint.gradient';
            obj.LHS = lhs;
%             obj.LHS = [obj.H,obj.constraint.gradient;obj.constraint.gradient',zeros(obj.m,obj.m)];
%            obj.LHS = [obj.H,obj.constraint.gradient';obj.constraint.gradient,zeros(obj.m,obj.m)];
        end
    end
end