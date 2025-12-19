classdef VolumeFunctionalRadius < handle

    properties (Access = private)
        mesh
        base
        test
    end

    properties (Access = private)
        baseFun
        totalVolume
    end

    methods (Access = public)
        function obj = VolumeFunctionalRadius(cParams)
            obj.init(cParams);
%             obj.createBaseFunction();
            obj.createTotalVolume();
        end

        function [J,dJ] = computeFunctionAndGradient(obj,x)
            xD = x.obtainDomainFunction();
%             xD = {x};
            J  = obj.computeFunction(xD{1});
            dJ{1} = obj.computeGradient(xD{1});
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.mesh = cParams.mesh;
%             obj.base = cParams.uMesh;
            obj.test = cParams.test;
            
            obj.baseFun = ConstantFunction.create(1,obj.mesh);
        end
        
%         function createBaseFunction(obj)
%             s.trial     = LagrangianFunction.create(obj.mesh,1,obj.test.order);
%             s.mesh      = obj.mesh;
%             riszFilter  = FilterLump(s);
%             f           = CharacteristicFunction.create(obj.base);
%             obj.baseFun = riszFilter.compute(f,2);
%         end

        function createTotalVolume(obj)           
            dV = obj.baseFun;
            V  = Integrator.compute(dV,obj.mesh,2);
            obj.totalVolume = V;
        end

        function J = computeFunction(obj,x)
            dV = obj.baseFun-pi*x.*x;
%             volume = Integrator.compute(dV,obj.mesh,2);
            volume = obj.totalVolume- pi*x.fValues'*x.fValues;
            J      = volume/obj.totalVolume;
        end

        function dJ = computeGradient(obj,x)
            dJ = copy(x);
            fValues = -2*pi*dJ.fValues;
            dJ.setFValues(fValues./obj.totalVolume);
            dJ.setFValues(fValues./norm(fValues));
        end
    end

    methods (Static, Access = public)
        function title = getTitleToPlot()
            title = 'Volume'; % Maybe a property in the future?
        end
    end
end