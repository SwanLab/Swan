classdef FilterPDEChambollePock < handle

    properties (Access = private)
        epsilon
        thetaRel
        tauG
        tauF
        mesh
        proxF
        proxG
        filterLump
        rho0
        sigma0
    end

    properties (Access = private)

    end

    methods (Access = public)

        function obj = FilterPDEChambollePock(cParams)
            obj.init(cParams);
            obj.createFilter();
            obj.createProximals();
            obj.rho0   = LagrangianFunction.create(obj.mesh, 1, 'P1');
            obj.sigma0 = LagrangianFunction.create(obj.mesh, 2, 'P1');
        end

        function rho = compute(obj,chi,quadType)
            proxG = @(rho) obj.proxG(rho,chi);
            [rho,sigma] = obj.solveWithChambollePockAlgorithm(obj.rho0,obj.sigma0,obj.proxF,proxG,obj.tauF,obj.tauG,obj.thetaRel);            
            obj.rho0 = rho;
            obj.sigma0 = sigma;
        end        

        function obj = updateEpsilon(obj,epsilon)
            if obj.hasEpsilonChanged(epsilon)
                obj.epsilon = epsilon;
            end 
        end

    end

    methods (Access = private)

        function init(obj,cParams)
           obj.mesh = cParams.mesh;
           obj.tauF = 0.0025;
           obj.tauG = 0.25;
           obj.thetaRel = 1;
           obj.epsilon  = 10;
        end

        function createFilter(obj)
            s.trial = LagrangianFunction.create(obj.mesh,1,'P1');
            s.mesh  = obj.mesh;
            obj.filterLump = FilterLump(s);
        end        

        function  createProximals(obj)
            %proxF = @(z)  proximalDroplet(z,tauF,k,alpha,ep);
            epsF    = 1; 
            obj.proxF = @(z)  obj.proximalEllipse(z,obj.tauF,epsF,obj.mesh);           
            tauGEps   = obj.tauG/obj.epsilon^2;
            obj.proxG = @(rho,chi) obj.proximalL2Projection(rho,chi,tauGEps);
        end        

        function [u,z] = solveWithChambollePockAlgorithm(obj,u0,z0,proxF,proxG,tauF,tauG,thetaRel)
            u  = u0;
            uN = u0;
            z   = z0;
            for kcp=1:10
                z      = proxF(z + tauF.*Grad(uN));
                uOld   = u;
                u      = proxG(u - tauG*Divergence(z));
                uN     = project(u + thetaRel.*(u - uOld),'P1');
            end
        end
        

      function s = proximalEllipse(obj,z,tau,alpha,m)
            A = [1 0; 0 1];
            I = eye(2);
            r = alpha^2/tau;
            dm = obj.createTensorFunction(A+r*I);
            s = LagrangianFunction.create(obj.mesh,2,'P1');
            M  = IntegrateLHS(@(u,v) DP(v,DP(dm,u)),s,s,m,'Domain');
            %LHS = diag(sum(M));
            LHS = M;
            F  = IntegrateRHS(@(v) DP(v,z),s,m,'Domain');
            sV = full(LHS\F);
            sV = reshape(sV,[s.ndimf,m.nnodes])';
            s.setFValues(sV);

        end

        function s = proximalDroplet(obj,z,tau,k,alpha,ep)
            s1 = z/(1+tau);
            zk  = z*k(:);
            z2 = sum(z.^2,2);
            tA   = tau/(alpha^2*(1+tau));
            tB   = sqrt((alpha^2-1)./(z2-zk.^2+ep));
            deltaS = tA*(z + (alpha^2-2)*zk.*k.' - tB.*((z2-2*zk.^2).*k.' + zk.*z));
            s      = s1 + (alpha*zk - sqrt(z2) > 0).*deltaS;
        end

      function Af = createTensorFunction(obj,A)
            op = @(xV) obj.evaluate(xV,A);
            Af = DomainFunction.create(op,obj.mesh,2);
        end

        function fV = evaluate(obj,xV,A)
            nGauss = size(xV,2);
            nElem  = obj.mesh.nelem;
            fV = repmat(A,[1 1 nGauss nElem]);
        end
        
        function rhoN = proximalL2Projection(obj,rho,Chi,tauG)
            rhoF = (1./(1+tauG)).*obj.filterLump.compute(rho,2);
            rhoChi = (1./(1+tauG)).*obj.filterLump.compute(Chi.*tauG,2);
            rhoN = project(rhoF + rhoChi,'P1');
        end        

        function itHas = hasEpsilonChanged(obj,eps)
            var   = abs(eps - obj.epsilon)/eps;
            itHas = var > 1e-15;
        end        


    end

end