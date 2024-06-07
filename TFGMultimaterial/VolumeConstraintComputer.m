classdef VolumeConstraintComputer < handle
    
    properties (Access = public)
        J
        dJ
    end

    properties (Access = private)
        area
        mesh
        vTar
        tfi
    end

    methods (Access = public)

        function obj = VolumeConstraintComputer(cParams)
            obj.init(cParams)
        end

        function [J,dJ] = computeFunctionAndGradient(obj,x)
            obj.computeCharacteristicFunction(x);
            J = obj.computeFunction();
            dJ = obj.computeGradient();
        end
    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh = cParams.mesh;
            obj.area  = cParams.area;
            obj.vTar = cParams.volumeTarget;
        end

        function computeCharacteristicFunction(obj,x)
            s.designVariable = x;
            s.m              = obj.mesh;
    
            charfun          = CharacteristicFunctionComputer(s); 
            [~,obj.tfi]         = charfun.computeFiandTfi();
        end

        function J = computeFunction(obj)
            surf = obj.area;
            charfunc = obj.tfi;

            V = surf*charfunc';
            J = (V/obj.vTar)-1;
        end

        function dJ = computeGradient(obj)
            dJ = 1./obj.vTar;
        end
    end

    methods (Static, Access = public)
        function title = getTitleToPlot()
            title = 'Volume';
        end
    end
end