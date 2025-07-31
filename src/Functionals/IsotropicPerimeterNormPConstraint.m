classdef IsotropicPerimeterNormPConstraint < handle

    properties (Access = private)
        mesh
        perTarget
        perimeter
    end
    
    methods (Access = public)
        function obj = IsotropicPerimeterNormPConstraint(cParams)
            obj.init(cParams);
        end
        
        function [J,dJ] = computeFunctionAndGradient(obj,x)
            [IP,dIP] = obj.perimeter.computeFunctionAndGradient(x);
            J        = obj.computeFunction(IP);
            dJ{1}    = obj.computeGradient(dIP{1});
        end  
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.mesh      = cParams.mesh;
            obj.perTarget = cParams.perTarget;
            obj.perimeter = IsotropicPerimeterNormPFunctional(cParams);
        end

        function J = computeFunction(obj,IP)
            pTar = obj.perTarget;
            J    = IP/pTar-1;
        end

        function dJ = computeGradient(obj,dIP)
            pTar = obj.perTarget;
            dJ   = dIP./pTar;
        end
    end

    methods (Static, Access = public)
        function title = getTitleToPlot()
            title = 'Iso Per p-norm constr';
        end
    end
end