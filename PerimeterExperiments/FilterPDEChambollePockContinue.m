classdef FilterPDEChambollePockContinue < handle

    properties (Access = private)
        epsilon
        thetaRel
        tauG
        tauF
        alpha
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

        function obj = FilterPDEChambollePockContinue(cParams)
            obj.init(cParams);
            obj.createFilter();
            obj.createProximals();
            obj.rho0   = LagrangianFunction.create(obj.mesh, 1, 'P1');
            obj.sigma0 = LagrangianFunction.create(obj.mesh, 2, 'P0');
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
           obj.thetaRel = 1;
           h = obj.mesh.computeMeanCellSize();
           obj.epsilon  = (4*h)^2;
           obj.tauF = 0.1*h; obj.epsilon;
           obj.tauG = 0.1*h; obj.epsilon;  
           obj.alpha = 0;%obj.epsilon;
           u = 45;
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

            %s.filterType = 'PDE';
            %s.mesh       = obj.mesh;
            %s.trial      = LagrangianFunction.create(obj.mesh,1,'P1');
            %f            = Filter.create(s);
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
            epsF    = obj.epsilon;%;1; 
            obj.proxF = @(z)  obj.proximalEllipse(z,obj.tauF,epsF,obj.mesh);           
            tauGEps   = obj.tauG;%/obj.epsilon^2;
            obj.proxG = @(rho,chi) obj.proximalL2Distance(rho,chi,tauGEps);
        end        

        function [u,z] = solveWithChambollePockAlgorithm(obj,u0,z0,proxF,proxG,tauF,tauG,thetaRel)
            u  = u0;
            uN = u0;
            z   = z0;
            err = 1; 
            iter = 1;
            while err > 1e-8
                z      = proxF(z + tauF.*Grad(uN));
                uOld   = u;
                divZ   = obj.computeDivergence(z);%Divergence(z)
                u      = proxG(project(u - tauG.*divZ,'P1'));
                uN     = project(u + thetaRel.*(u - uOld),'P1');
                err = Norm(u-obj.rhoExact,'L2')/Norm(obj.rhoExact,'L2')
                error(iter) = err;
                iter = iter+1;
            end
            plot(log10(error))
            plot(u)
            plot(obj.rhoExact)
        end

        function rhoN = computeDivergence(obj,z)
            a = obj.alpha;            
            m = obj.mesh;
            rhoN = LagrangianFunction.create(obj.mesh,1,'P1');
            M  = IntegrateLHS(@(u,v) DP(v,u),rhoN,rhoN,m,'Domain');
            K  = IntegrateLHS(@(u,v) (a^2).*DP(Grad(v),Grad(u)),rhoN,rhoN,m,'Domain');
            F = IntegrateRHS(@(v) DP(Grad(v),z),rhoN,m,'Domain');
            rhoV = full((K+M)\F);
            rhoV = reshape(rhoV,rhoN.ndimf,[])';
            rhoN.setFValues(rhoV);           
        end
        

      function s = proximalEllipse(obj,z,tau,eps,m)
            %A = (obj.CGlobal);%inv([1 0; 0 1]);
            I = eye(2);
            %r = eps^2/tau;
 %          dm = obj.createTensorFunction(A+r*I);
            r = tau/eps;
            dm = obj.createTensorFunction(I+r*I);
            s = LagrangianFunction.create(obj.mesh,2,'P0');
            M  = IntegrateLHS(@(u,v) DP(v,DP(dm,u)),s,s,m,'Domain');
           % LHS = diag(sum(M));
            LHS = M;
            F  = IntegrateRHS(@(v) DP(v,z),s,m,'Domain');
            sV = full(LHS\F);
            sV = reshape(sV,s.ndimf,[])';
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
        
        function rhoN = proximalL2Distance(obj,rho,chi,tauG)
            rhoN = LagrangianFunction.create(obj.mesh,1,'P1');
            m    = obj.mesh;            
            e = obj.epsilon;
            a = obj.alpha;
            M  = IntegrateLHS(@(u,v) (1+e/tauG).*DP(v,u),rhoN,rhoN,m,'Domain');
            K  = IntegrateLHS(@(u,v) (e*a^2/tauG).*DP(Grad(v),Grad(u)),rhoN,rhoN,m,'Domain');
            F1 = IntegrateRHS(@(v) DP(v,chi+e/tauG*rho),rhoN,m,'Domain');
            F2 = IntegrateRHS(@(v) e*a^2/tauG.*DP(Grad(v),Grad(rho)),rhoN,m,'Domain');
            rhoV = full((K+M)\(F1+F2));
            rhoV = reshape(rhoV,rhoN.ndimf,[])';
            rhoN.setFValues(rhoV);
        end        

        function itHas = hasEpsilonChanged(obj,eps)
            var   = abs(eps - obj.epsilon)/eps;
            itHas = var > 1e-15;
        end        


    end

end