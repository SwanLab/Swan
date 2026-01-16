classdef FilterPDEChambollePockContinue < handle

    properties (Access = private)
        epsilon
        thetaRel
        tauG
        tauF
        eta
        k
        mesh
        proxF
        proxG
        filterLump
        rho0
        sigma0
        rhoExact
        CGlobal
        alpha
        beta
        tol
        MassRho
        MassSigma
    end

    properties (Access = private)

    end

    methods (Access = public)

        function obj = FilterPDEChambollePockContinue(cParams)
            obj.init(cParams);
            %obj.createFilter();
            obj.createProximals();
            obj.rho0   = LagrangianFunction.create(obj.mesh, 1, 'P1');
            obj.sigma0 = LagrangianFunction.create(obj.mesh, 2, 'P0');
        end
    
        function rho = compute(obj,chi,quadType)
         %   obj.rhoExact = obj.computeExactSolution(chi);
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
           obj.epsilon  = 0.1;%(4*h)^2;
           obj.tauF = 0.1*h;% obj.epsilon; 
           obj.tauG = 0.1*h;% obj.epsilon;  
           obj.eta = 0;%obj.epsilon; %0
           obj.k = [1 1]';
           obj.k = obj.k/norm(obj.k);
           obj.alpha = 1;
           obj.beta  = obj.epsilon;0.001;
           obj.tol   = 0.001;%99;%1e0;
           u = 85;
           gamma = 90;
           CAnisotropic = [tand(u),0;0,1/tand(u)];
           R = [cosd(gamma),-sind(gamma)
               sind(gamma), cosd(gamma)];
           obj.CGlobal = R*CAnisotropic*R';
           rhoN = LagrangianFunction.create(obj.mesh,1,'P1');
           obj.MassRho = IntegrateLHS(@(u,v) DP(v,u),rhoN,rhoN,obj.mesh,'Domain');
        end

        function rho = computeExactSolution(obj,chi)
            % A       =  ConstantFunction.create(obj.CGlobal,obj.mesh);
            % s.A     = A;
            % s.trial = LagrangianFunction.create(obj.mesh,1,'P1');
            % s.mesh = obj.mesh;
            % s.filterType = 'PDE';
            % s.boundaryType = 'Neumann';
            % s.metric = 'Anisotropy';
            % f = Filter.create(s);


            s.alpha = obj.alpha;
            s.beta  = obj.beta;
            s.mesh  = obj.mesh;
            s.theta = obj.k;
            s.filterType = 'Segment';
            s.tol0  = obj.tol;
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
            %obj.proxF = @(z)  obj.proxEllipse(z,obj.tauF,epsF,obj.mesh);           
            obj.proxF = @(z)  obj.proxSegment(z,obj.tauF,epsF,obj.mesh);           
            tauGEps   = obj.tauG;%/obj.epsilon^2;
            obj.proxG = @(rho,chi) obj.proximalL2Distance(rho,chi,tauGEps);
        end          

        function [u,z] = solveWithChambollePockAlgorithm(obj,u0,z0,proxF,proxG,tauF,tauG,thetaRel)
            u  = u0;
            uN = u0;
            z   = z0;
            stopCrit = 1; 
            iter = 1;
            while stopCrit > obj.tol
                z      = proxF(z + tauF.*Grad(uN));
                uOld   = u;
                divZ   = obj.computeDivergence(z);%Divergence(z)

                uStar = copy(u);
                uStar.setFValues(u.fValues - tauG.*divZ.fValues);
                u      = proxG(uStar);
                
                uN.setFValues(u.fValues + thetaRel.*(u.fValues - uOld.fValues));
                stopCrit = Norm(u-uOld,'H1',obj.eta)/(Norm(uOld,'H1',obj.eta)+1e-5);
                error(iter) = stopCrit;
                iter = iter+1
            end
           % plot(log10(error))
           % plot(u)
           % plot(obj.rhoExact)
        end

        function rhoN = computeDivergence(obj,z)
            a = obj.eta;            
            m = obj.mesh;
           % rhoN = LagrangianFunction.create(obj.mesh,1,'P1');
          %  M  = IntegrateLHS(@(u,v) DP(v,u),rhoN,rhoN,m,'Domain');
            %K  = IntegrateLHS(@(u,v) (a^2).*DP(Grad(v),Grad(u)),rhoN,rhoN,m,'Domain');
            F = IntegrateRHS(@(v) DP(Grad(v),z),obj.rho0,m,'Domain');
            %LHS = M+K;
            LHS = obj.MassRho;
           % LHS = diag(sum(M));
            rhoN = copy(obj.rho0);
            rhoV = full(LHS\F);
            rhoV = reshape(rhoV,rhoN.ndimf,[])';
            rhoN.setFValues(rhoV);           
        end
        

      function s = proxEllipse(obj,z,tau,eps,m)
            Ainv = inv(obj.CGlobal);%inv([1 0; 0 1]);
            I = eye(2);
            %r = eps^2/tau;
 %          dm = obj.createTensorFunction(A+r*I);
            r = tau/eps;
            dm = obj.createTensorFunction(r*Ainv+I);
            s = LagrangianFunction.create(obj.mesh,2,'P0');
            M  = IntegrateLHS(@(u,v) DP(v,DP(dm,u)),s,s,m,'Domain');
           % LHS = diag(sum(M));
            LHS = M;
            F  = IntegrateRHS(@(v) DP(v,z),s,m,'Domain');
            sV = full(LHS\F);
            sV = reshape(sV,s.ndimf,[])';
            s.setFValues(sV);
      end

      function s = proxSegment(obj,z,tau,eps,m)
            a = obj.alpha;
            b = obj.beta;
            k = obj.k;
            %r = eps^2/tau;
 %          dm = obj.createTensorFunction(A+r*I);
            ra = 1/(tau/(eps*a^2)+1);
            rb = 1/(tau/(eps*b^2)+1);
            z = project(z,'P0');
            sV  = (ra*max(0,z.fValues*k) + rb*min(0,z.fValues*k))*k';


            s = LagrangianFunction.create(obj.mesh,2,'P0');
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
           % A       =  ConstantFunction.create(obj.CGlobal,obj.mesh);
            e = obj.epsilon;
            % a = obj.eta;
            rhoN = copy(obj.rho0);
            % m    = obj.mesh;            
           
            % M  = IntegrateLHS(@(u,v) (1+e/tauG).*DP(v,u),rhoN,rhoN,m,'Domain');
            % K  = IntegrateLHS(@(u,v) (e*a^2/tauG).*DP(Grad(v),Grad(u)),rhoN,rhoN,m,'Domain');
            % F1 = IntegrateRHS(@(v) DP(v,chi+e/tauG*rho),rhoN,m,'Domain');
            % F2 = IntegrateRHS(@(v) e*a^2/tauG.*DP(Grad(v),Grad(rho)),rhoN,m,'Domain');
            % rhoV = full((K+M)\(F1+F2));
            

            rhoV = (chi.fValues+e/tauG*rho.fValues)/(1+e/tauG);
            rhoV = reshape(rhoV,rhoN.ndimf,[])';
            rhoN.setFValues(rhoV);

        end        

        function itHas = hasEpsilonChanged(obj,eps)
            var   = abs(eps - obj.epsilon)/eps;
            itHas = var > 1e-15;
        end        


    end

end