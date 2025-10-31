classdef FilterPDEChambollePock < handle

    properties (Access = public)

    end

    properties (Access = private)

    end

    properties (Access = private)

    end

    methods (Access = public)

        function obj = FilterPDEChambollePock(cParams)
            obj.init(cParams)

        end

        function xF = compute(obj,fun,quadType)
            xF = LagrangianFunction.create(obj.mesh, fun.ndimf, obj.trial.order);
         end        

    end

    methods (Access = private)

        function init(obj)

        end

        function

        end

    end

end