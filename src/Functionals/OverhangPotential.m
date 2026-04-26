classdef OverhangPotential < handle

    properties (Access = private)
        mesh
        charFun
        epsilon
        vec
        th
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
            k     = obj.vec;
            gRhoN = sqrt(DP(Grad(rho),Grad(rho)));
            gRhoU = Grad(rho)./(gRhoN + 1e-6);
            dirDer = DP(gRhoU,k);
            J1  = obj.computeMinimumSquaresTerm(rho);
            J2  = obj.computeRegularizationTerm(gRhoN,dirDer);

            rhs1 = obj.computeRHSMinimumSquares(rho);
            rhs2 = obj.computeRHSRegTerm1(rho,gRhoN,gRhoU,dirDer);
            rhs3 = obj.computeRHSRegTerm2(rho,gRhoN,gRhoU,dirDer);

            J = J1 + J2;
            rhs = rhs1 + rhs2 + rhs3;
            if isempty(obj.value0)
                obj.value0 = J;
            end
            J = J/obj.value0;
            dJ{1} = copy(rho);
%             dJ{1}.setFValues(rhs./(obj.value0.*obj.M));

            h = obj.mesh.computeMeanCellSize();
            f = @(v,u) v.*u;
            g = @(v,u) DP(Grad(v),Grad(u));
            M = IntegrateLHS(f,rho,rho,obj.mesh,'Domain',2);
            K = IntegrateLHS(g,rho,rho,obj.mesh,'Domain',2);
            LHS = M+h^2*K;
            dJ{1}.setFValues(LHS\rhs);
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.mesh    = cParams.mesh;
            obj.charFun = cParams.chi;
            obj.epsilon = cParams.epsilon;
            obj.vec     = cParams.k;
            obj.th      = cParams.theta;
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

        function J = computeRegularizationTerm(obj,gRhoN,dirDer)
            theta  = obj.th;
            e      = obj.epsilon;
            f      = gRhoN./(dirDer - cos(theta) + 1e-6);
            int    = (f.^2).*(max(dirDer-cos(theta),0)).^2;
            J      = (e^2/2)*Integrator.compute(int,obj.mesh,3);
        end

        function rhs = computeRHSMinimumSquares(obj,rho)
            chi  = obj.charFun;
            rhs1 = IntegrateRHS(@(v) rho.*v,rho,obj.mesh,'Domain',2);
            rhs2 = -IntegrateRHS(@(v) chi.*v,rho,obj.mesh,'Domain',2);
            rhs  = rhs1 + rhs2;
        end

        function rhs = computeRHSRegTerm1(obj,rho,gRhoN,gRhoU,dirDer)
            e = obj.epsilon;
            k = obj.vec;
            theta = obj.th;
            constr = dirDer - cos(theta) + 1e-6;
            maxF = max(dirDer - cos(theta),0).^2;
            gV    = @(v) Grad(v);
            gradV = @(v) DomainFunction.create(@(xV) squeezeParticular(gV(v).evaluate(xV),1),obj.mesh,gRhoN.ndimf);
            num1 = @(v) constr.*DP(Grad(rho),gradV(v));
            num2 = @(v) gRhoN.*DP(k,gradV(v)-Grad(rho).*DP(gRhoU,gradV(v)./(gRhoN + 1e-6)));
            int  = @(v) maxF.*(num1(v)-num2(v))./(constr.^3);
            rhs  = e^2.*IntegrateRHS(int,rho,obj.mesh,'Domain',3);
        end

        function rhs = computeRHSRegTerm2(obj,rho,gRhoN,gRhoU,dirDer)
            e = obj.epsilon;
            k = obj.vec;
            theta = obj.th;
            constr = dirDer - cos(theta) + 1e-6;
            maxF = max(dirDer - cos(theta),0);
            gV    = @(v) Grad(v);
            gradV = @(v) DomainFunction.create(@(xV) squeezeParticular(gV(v).evaluate(xV),1),obj.mesh,gRhoN.ndimf);
            f1   = (DP(Grad(rho),Grad(rho)))./(constr.^2);
            f2   = @(v) DP(k,gradV(v)./(gRhoN + 1e-6) - gRhoU.*DP(gRhoU,gradV(v)./(gRhoN + 1e-6)));
            int = @(v) maxF.*f1.*f2(v);
            rhs  = e^2.*IntegrateRHS(int,rho,obj.mesh,'Domain',3);
        end
    end

    methods (Static, Access = public)
        function title = getTitleToPlot()
            title = 'Ov potential';
        end
    end
end