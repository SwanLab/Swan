classdef VolumeFunctionalParameter < handle

    properties (Access = private)
        mesh
        base
        test
    end

    properties (Access = private)
        baseFun
        totalVolume
        volume
        gradJ
    end

    methods (Access = public)
        function obj = VolumeFunctionalParameter(cParams)
            obj.init(cParams);
%             obj.createBaseFunction();
            obj.createTotalVolume();
        end

        function [J,dJ] = computeFunctionAndGradient(obj,x)
            xD = x.obtainDomainFunction();
%             xD = {x};
            J  = obj.computeFunction(xD{1});
            dJ = obj.computeGradient(xD{1});
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.mesh = cParams.mesh;
%             obj.base = cParams.uMesh;
            obj.test = cParams.test;            
            obj.baseFun = ConstantFunction.create(1,obj.mesh);
            % switch cParams.geomType
            %     case 'Circle'
            %         obj.volumeHole = @(x) pi*x.fValues'*x.fValues;
            %         obj.gradJ      = @(x) -2*pi*x.fValues;
            %     case 'Square'
            %         obj.volumeHole = @(x) 4*x.fValues'*x.fValues;
            %         obj.gradJ      = @(x) -8*x.fValues;
            %     case 'Lattice'
            %         obj.volumeHole = @(x) 4*x.fValues'*x.fValues;
            %         obj.gradJ      = @(x) -8*x.fValues;
            % end
            obj.volume = cParams.volume;
            obj.gradJ  = cParams.gradJ;
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
            volume = obj.volume(x);
            J      = volume/obj.totalVolume;
           % rho = (L^2-x)./L^2;
%            J   = Integrator.compute(rho,obj.mesh,2)/obj.totalVolume;
        end

        function dJ = computeGradient(obj,x)
            % dJ = copy(x);
            s.mesh  = x.mesh;
            s.order = x.order;
            nparam  = x.ndimf; 
            fValues = obj.gradJ(x);
            for i = 1:nparam
                %s.fValues = fValues(:,i)/norm(fValues(:,i));
                 s.fValues = fValues(:,i)/obj.totalVolume;
                dJ{i} = LagrangianFunction(s);
            end

%             fValues = -8*dJ.fValues;
            % fValues = obj.gradJ(x);
%             dJ.setFValues(fValues./(obj.totalVolume));
%              dJ.setFValues(fValues);
            % dJ.setFValues(fValues./norm(fValues));
        end
    end

    methods (Static, Access = public)
        function title = getTitleToPlot()
            title = 'Volume'; % Maybe a property in the future?
        end
    end
end
