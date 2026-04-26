classdef OverhangPotential < handle

    properties (Access = private)
        mesh
        charFun
        M
    end

    properties (Access = private)
        value0
    end

    methods (Access = public)
        function obj = OverhangPotential(cParams)
            obj.init(cParams);
        end

        function [J,dJ] = computeFunctionAndGradient(obj,x)
            xD  = x.obtainDomainFunction();
            rho = xD{1};
            J1  = obj.computeMinimumSquaresTerm(rho);

            rhs1 = obj.computeRHSMinimumSquares(rho);

            J = J1;
            rhs = rhs1;
            if isempty(obj.value0)
                obj.value0 = J;
            end
            J = J/obj.value0;
            dJ{1} = copy(rho);
            dJ{1}.setFValues(rhs./(obj.value0.*obj.M));
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.mesh    = cParams.mesh;
            obj.charFun = cParams.chi;
            obj.M       = obj.createLumpedMass(cParams);
        end

        function M = createLumpedMass(obj,cParams)
            f = @(v,u) v.*u;
            M = IntegrateLHS(f,cParams.trial,cParams.trial,obj.mesh,'Domain',2);
            M = sum(M,2);
        end

        function J = computeMinimumSquaresTerm(obj,rho)
            chi = obj.charFun;
            int1 = Integrator.compute(rho.*rho,obj.mesh,2);
            int2 = -Integrator.compute(chi.*(rho.*2),obj.mesh,2);
            int3 = Integrator.compute(chi.*chi,obj.mesh,2);
            J    = 0.5*(int1+int2+int3);
        end

        function rhs = computeRHSMinimumSquares(obj,rho)
            chi  = obj.charFun;
            rhs1 = IntegrateRHS(@(v) rho.*v,rho,obj.mesh,'Domain',2);
            rhs2 = -IntegrateRHS(@(v) chi.*v,rho,obj.mesh,'Domain',2);
            rhs  = rhs1 + rhs2;
        end
    end

    methods (Static, Access = public)
        function title = getTitleToPlot()
            title = 'Ov potential';
        end
    end
end