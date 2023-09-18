classdef IPMSymmetricLHSComputer < handle

    properties (Access = public)
        LHS
    end

    properties (Access = private)
        H
        constraint
        m
    end

    properties (Access = private)
        nX
        nSlack
        hessian
        upperZ
        lowerZ
        diagonaldL
        diagonaldU
    end

    methods (Access = public)
        function obj = IPMSymmetricLHSComputer(cParams)
            obj.init(cParams);
        end
        function LHS = compute(obj)
            LHS = obj.computeLHS();
        end
    end

    methods (Access = private)

        function init(obj,cParams)
            obj.H = cParams.H;
            obj.constraint = cParams.constraint;
            obj.m = cParams.m;
        end

        function lhs = computeLHS(obj)
            sizeH = size(obj.H,1);
            sizeg = obj.m;
            lhs = zeros(sizeH+sizeg);
            lhs(1:sizeH,1:sizeH) = obj.H;
            lhs(1:sizeH,sizeH+1:sizeH+sizeg) = obj.constraint.gradient;
            lhs(sizeH+1:sizeH+sizeg,1:sizeH) = obj.constraint.gradient';
        end
    end
end