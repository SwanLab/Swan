classdef OverhangPotential < handle

    properties (Access = private)
        mesh
        charFun
        epsilon
        vec
        th
        LHS
        value0
        tol
    end

    methods (Access = public)
        function obj = OverhangPotential(cParams)
            obj.init(cParams);
        end

        function [J,dJ] = computeFunctionAndGradient(obj,x)
            xD      = x.obtainDomainFunction();
            rho     = xD{1};
            k       = obj.vec;
            gRho    = Grad(rho);
            gRhoN   = sqrt(DP(Grad(rho),Grad(rho)));
            dirDer  = DP(gRho,k);
            theta   = obj.th;
            constr  = dirDer - gRhoN.*cos(theta);
            obj.tol = 1e-5;
            if isempty(obj.tol)
                obj.tol = 1e-30;
            end
            J1      = obj.computeMinimumSquaresTerm(rho);
            J2      = obj.computeRegularizationTerm(gRho,constr);

            rhs1 = obj.computeRHSMinimumSquares(rho);
            rhs2 = obj.computeRHSRegTerm1(rho,gRho,constr);
            rhs3 = obj.computeRHSRegTerm2(rho,gRho,gRhoN,constr);
            rhs4 = obj.computeRHSRegTerm3(rho,gRho,gRhoN,constr);

            J = J1 + J2;
            rhs = rhs1 + rhs2 + rhs3 + rhs4;
            J = J/obj.value0;
            dJ{1} = copy(rho);
            dJ{1}.setFValues(obj.LHS\(rhs./(obj.value0)));

        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.mesh    = cParams.mesh;
            obj.charFun = cParams.chi;
            obj.epsilon = cParams.epsilon;
            obj.vec     = cParams.k;
            obj.th      = cParams.theta;
            obj.LHS     = obj.createLHS(cParams);
            obj.value0  = obj.computeMinimumSquaresTerm(obj.charFun+obj.charFun);
        end

        function MK = createLHS(obj,cParams)
            h = obj.mesh.computeMeanCellSize();
            f = @(v,u) v.*u;
            M = IntegrateLHS(f,cParams.trial,cParams.trial,obj.mesh,'Domain',2);
            g = @(v,u) DP(Grad(v),Grad(u));
            K = IntegrateLHS(g,cParams.trial,cParams.trial,obj.mesh,'Domain',2);
            MK = M + h^2.*K;
        end

        function J = computeMinimumSquaresTerm(obj,rho)
            chi = obj.charFun;
            int1 = Integrator.compute(rho.*rho,obj.mesh,2);
            int2 = -2*Integrator.compute(chi.*rho,obj.mesh,2);
            int3 = Integrator.compute(chi.*chi,obj.mesh,2);
            J    = 0.5*(int1+int2+int3);
        end

        function J = computeRegularizationTerm(obj,gRho,constr)
            e      = obj.epsilon;
            f      = DP(gRho,gRho)./(constr + obj.tol).^2;
            int    = f.*(max(constr,0)).^2;
            J      = (e^2/2)*Integrator.compute(int,obj.mesh,3);
        end

        function rhs = computeRHSMinimumSquares(obj,rho)
            chi  = obj.charFun;
            rhs1 = IntegrateRHS(@(v) rho.*v,rho,obj.mesh,'Domain',2);
            rhs2 = -IntegrateRHS(@(v) chi.*v,rho,obj.mesh,'Domain',2);
            rhs  = rhs1 + rhs2;
        end

        function rhs = computeRHSRegTerm1(obj,rho,gRho,constr)
            e    = obj.epsilon;
            maxF = max(constr,0).^2;
            int  = @(v) DP(gRho,Grad(v)).*(maxF./(constr + obj.tol).^2);
            rhs  = e^2.*IntegrateRHS(int,rho,obj.mesh,'Domain',3);
        end

        function rhs = computeRHSRegTerm2(obj,rho,gRho,gRhoN,constr)
            e = obj.epsilon;
            k = obj.vec;
            theta = obj.th;
            maxF = max(constr,0);
            den = (constr + obj.tol).^2;
            funV = @(v) DP(gRhoN.*k,Grad(v)) - DP(cos(theta).*gRho,Grad(v));
            int = @(v) gRhoN.*(maxF./den).*funV(v);
            rhs  = e^2.*IntegrateRHS(int,rho,obj.mesh,'Domain',3);
        end

        function rhs = computeRHSRegTerm3(obj,rho,gRho,gRhoN,constr)
            e = obj.epsilon;
            k = obj.vec;
            theta = obj.th;
            maxF = max(constr,0).^2;
            den = (constr + obj.tol).^3;
            funV = @(v) DP(gRhoN.*k,Grad(v)) - DP(cos(theta).*gRho,Grad(v));
            int = @(v) gRhoN.*(maxF./den).*funV(v);
            rhs  = -e^2.*IntegrateRHS(int,rho,obj.mesh,'Domain',3);
        end
    end

    methods (Static, Access = public)
        function title = getTitleToPlot()
            title = 'Ov potential';
        end
    end
end