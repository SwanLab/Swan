classdef TopologicalDerivativeComputer < handle

    properties (Access = public)
        dt
        dC
        strain
        TD
        tgamma
    end

    properties (Access = private)
        tensor
        %strain
        meshSeba
        p
        t
        mesh
        psi
        designVariable
        penalization
        penalty
        volfrac
        auglag
        max_vol
        energy0
        voltarget
        nMat
        mat
        gamma
        
        tE
        beta
        alpha
        la0 
        mu0
        young
        stress
        U
        volume
        swanStress
    end

    methods (Access = public)

        function obj = TopologicalDerivativeComputer(cParams)
            obj.init(cParams)
            obj.computeGamma();
            obj.computeTgamma();
            obj.computeBetaAndAlpha();
            obj.computeLameParameters();
            obj.computeStressAndStrain();
            obj.computeTopologicalDerivativeAndDC();
        end
    end

    methods (Access = private)
        
        function init(obj,cParams)
            obj.p = cParams.mesh.coord';
            obj.t = cParams.mesh.connec';
            % obj.t(4,:) = 1;
            obj.mesh = cParams.mesh;
            obj.nMat = cParams.nMat;
            obj.mat{1} = cParams.mat.A;
            obj.mat{2} = cParams.mat.B;
            obj.mat{3} = cParams.mat.C;
            obj.mat{4} = cParams.mat.D;
            %obj.psi = cParams.psi;   
            obj.designVariable = cParams.designVariable;
            obj.U = cParams.U;
            obj.volume = cParams.volume;
        
        end

        function computeGamma(obj)
            E1 = obj.mat{1}.young;
            E2 = obj.mat{2}.young;
            E3 = obj.mat{3}.young;
            E4 = obj.mat{4}.young;

            obj.young = [E1 E2 E3 E4];

            obj.gamma(1) = E1/E1;
            obj.gamma(2) = E2/E1;
            obj.gamma(3) = E3/E1;
            obj.gamma(4) = E4/E1;    
        end

        function computeTgamma(obj)
            %s.psi = obj.psi;
            s.designVariable = obj.designVariable;
            s.m = obj.mesh;
            
            charfun = CharacteristicFunctionComputer(s); 
            [~,tfi] = charfun.computeFiandTfi();

            obj.tgamma = obj.gamma*tfi; %Mixed formulation method
            E1 = obj.mat{1}.young;
            obj.tE = E1*obj.tgamma; 
        end

        function computeBetaAndAlpha(obj)
             
            nu = obj.mat{1}.nu;
            obj.beta = (1+nu)/(1-nu); 
            obj.alpha = (3-nu)/(1+nu);
        end

        function computeLameParameters(obj)
            nu = obj.mat{1}.nu;
            E1 = obj.mat{1}.young;
        
            obj.la0 = nu.*E1./((1+nu).*(1-2.*nu)); obj.mu0 = E1./(2.*(1+nu)); % plane strain
            obj.la0 = 2.*obj.mu0.*obj.la0./(obj.la0+2.*obj.mu0); % plane stress 

        end

        function computeStressAndStrain(obj)
            [ux,uy]=pdegrad(obj.p,obj.t,obj.U); % solution gradient
            obj.strain=[ux(1,:);(ux(2,:)+uy(1,:))/2;uy(2,:)]; % strain
            e = obj.strain;
            obj.tensor = [obj.la0+2*obj.mu0 0 obj.la0; 0 2*obj.mu0 0; obj.la0 0 obj.la0+2*obj.mu0];
            sSwan = obj.tensor*e;
            % effective stress
            tgamma3 = [obj.tgamma;obj.tgamma;obj.tgamma]; 
            obj.swanStress = sSwan.*tgamma3; 
        end

        function computeTopologicalDerivativeAndDC(obj)
            sSwan = obj.swanStress;
            C = obj.tensor;
            e = obj.strain;
            E = obj.young;
            a = obj.alpha;
            b = obj.beta;

            for i = 1:obj.nMat
                for j = 1:obj.nMat
                    g  = E(j)/E(i);
                    c1 = 0.5*((1-g)./(1+a*g))./obj.tE;
                    c2 = c1.*((g.*(a-2*b)-1)./(1+b*g));

                    for k=1:size(sSwan,2)
                        coefMatrix = [4*c1(k)+c2(k) 0 c2(k); 0 8*c1(k) 0; c2(k) 0 4*c1(k)+c2(k)]; 
                        obj.dC(:,:,k,i,j) = C*coefMatrix*C;
                    end
                end
            end

        end
        
    end
        
end