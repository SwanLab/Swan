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
        rhoExact
        CGlobal
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
            obj.rhoExact = obj.computeExactSolution(chi);
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
           obj.tauF = 2.5*1e4;
           obj.tauG = 5*1e-5;
           obj.thetaRel = 2;
           obj.epsilon  = (4*obj.mesh.computeMeanCellSize())^2;
           u = 85;
           alpha = 90;
           CAnisotropic = [tand(u),0;0,1/tand(u)];
           R = [cosd(alpha),-sind(alpha)
               sind(alpha), cosd(alpha)];
           obj.CGlobal = R*CAnisotropic*R';
        end

        function rho = computeExactSolution(obj,chi)


            A       =  ConstantFunction.create(obj.CGlobal,obj.mesh);
            s.A     = A;
            s.trial = LagrangianFunction.create(obj.mesh,1,'P1');
            s.mesh = obj.mesh;
            s.filterType = 'PDE';
            s.boundaryType = 'Neumann';
            s.metric = 'Anisotropy';
            f = Filter.create(s);

           s.filterType = 'PDE';
            s.mesh       = obj.mesh;
            s.trial      = LagrangianFunction.create(obj.mesh,1,'P1');
            f            = Filter.create(s);
            f.updateEpsilon(obj.epsilon);
            rho = f.compute(chi,2);        
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
            for iter = 1:5000
                z      = proxF(z + tauF.*Grad(uN));
                uOld   = u;
                u      = proxG(u - tauG*Divergence(z));
                uN     = project(u + thetaRel.*(u - uOld),'P1');
                err = Norm(u-obj.rhoExact,'L2')/Norm(obj.rhoExact,'L2')
                error(iter) = err;
            end
            plot(error)
            plot(u)
            plot(obj.rhoExact)
        end
        

      function s = proximalEllipse(obj,z,tau,alpha,m)
            A = (obj.CGlobal);%inv([1 0; 0 1]);
            I = eye(2);
            r = alpha^2/tau;
            dm = obj.createTensorFunction(A+r*I);
            s = LagrangianFunction.create(obj.mesh,2,'P1');
            M  = IntegrateLHS(@(u,v) DP(v,DP(dm,u)),s,s,m,'Domain');
            LHS = diag(sum(M));
            %LHS = M;
            F  = IntegrateRHS(@(v) r*DP(v,z),s,m,'Domain');
            sV = full(LHS\F);
            sV = reshape(sV,[s.ndimf,m.nnodes])';
            s.setFValues(sV);

        end

        % function s = proximalDroplet(obj,z,tau,k,alpha,ep)
        %     s1 = z/(1+tau);
        %     zk  = z*k(:);
        %     z2 = sum(z.^2,2);
        %     tA   = tau/(alpha^2*(1+tau));
        %     tB   = sqrt((alpha^2-1)./(z2-zk.^2+ep));
        %     deltaS = tA*(z + (alpha^2-2)*zk.*k.' - tB.*((z2-2*zk.^2).*k.' + zk.*z));
        %     s      = s1 + (alpha*zk - sqrt(z2) > 0).*deltaS;
        % end

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