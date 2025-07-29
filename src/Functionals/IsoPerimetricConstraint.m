classdef IsoPerimetricConstraint < handle

    properties (Access = private)
        mesh
        C
        isoPer
    end
    
    methods (Access = public)
        function obj = IsoPerimetricConstraint(cParams)
            obj.init(cParams);
        end
        
        function [J,dJ] = computeFunctionAndGradient(obj,x)
            [IP,dIP] = obj.isoPer.computeFunctionAndGradient(x);
            J        = obj.computeFunction(IP);
            dJ{1}    = obj.computeGradient(dIP{1});
        end  
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.mesh   = cParams.mesh;
            obj.C      = cParams.C;
            obj.isoPer = IsoPerimetricFunctional(cParams);
        end

        function J = computeFunction(obj,IP)
            CTar = obj.C;
            J    = IP/CTar-1;
        end

        function dJ = computeGradient(obj,dIP)
            CTar = obj.C;
            dJ   = dIP./CTar;
        end
    end

    methods (Static, Access = public)
        function title = getTitleToPlot()
            title = 'Isoperimetric constraint';
        end
    end
end

