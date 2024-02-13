classdef HomogenizedVarComputerFromVademecum ...
        < HomogenizedVarComputer
    
    properties (Access = public)
        Ctensor
        density
        Ptensor
        rotator
    end
    
    properties (Access = private)
       vademecum
       superEllipseExponent
       PpRef
       dPpRef
       P2Ref
       dP2Ref
       superEllipseExp
       xiSymmetry
       pNorm
       xV
    end
    
    properties (Access = private)
       fileName 
    end
    
    methods (Access = public)
        
        function obj = HomogenizedVarComputerFromVademecum(cParams)
            obj.init(cParams);
            obj.loadVademecum();
            obj.loadHomogenizedVariables();
            obj.createRotator();
        end
        

        function computePtensor(obj,x,pNorm)
            obj.pNorm = pNorm;
            obj.obtainReferenceAmplificatorTensor(x,pNorm);
            alpha = x{3};
            obj.rotateAmplificatorTensor(alpha);
        end
        
        function computeDensity(obj,x)
            mx = x{1};
            my = x{2};
            [rho,drho] = obj.density.compute([mx,my]);
            obj.rho = rho;
            obj.drho{1} = drho(:,:,1);
            obj.drho{2} = drho(:,:,2);
        end
        
        function fP = addPrintableVariables(obj,x)
            obj.superEllipseExp = obj.computeSuperEllipseExponent(x);
            obj.computeXiSymmetry(x);
            fP{1}.value = obj.rho;
            fP{2}.value = obj.superEllipseExp;
            fP{3}.value = obj.xiSymmetry;
            pMax = 32;
            fP{4}.value = obj.computePcomponents(pMax,x);
        end
        
        function PpV = computePcomponents(obj,p,x)
           % xR = obj.reshapeDesignVariable(x);
           % mx = xR{1};
          %  my = xR{2};
            %[~,~,Pp,~] = obj.Ptensor.compute([mx,my],p);
            Pp = obj.PpRef;
            Ppv = (abs(Pp)).^(1/(p/2))';
            [PpV,ind] = max(Ppv');
            PpV = PpV';
        end
        
        function fP = createPrintVariables(obj)
            fP{1}.type = 'ScalarGauss';
            fP{2}.type = 'ScalarNodal';
            fP{3}.type = 'ScalarNodal';
            fP{4}.type = 'ScalarGauss';
            fP{1}.name = 'DensityGauss';
            fP{2}.name = 'SuperEllipseExponent';
            fP{3}.name = 'XiSymmetry';
            fP{4}.name = 'Pmax';
        end
        
    end
    
    methods (Access = private)
        

        function loadHomogenizedVariables(obj)
            obj.Ctensor = obj.vademecum.Ctensor;
            obj.Ptensor = obj.vademecum.Ptensor;
            obj.density = obj.vademecum.density;
        end
        
        function createSuperEllipseExponent(obj,m1,m2)
            switch obj.fileName
                case 'SuperEllipseQMax'
                    s.type = 'Given';
                    s.q    = 32*ones(size(m1));
                case 'SuperEllipseQ2' 
                    s.type = 'Given';
                    s.q    = 2*ones(size(m1));
                case 'SuperEllipseQOptAnalytic'
                    s.type = 'Optimal';
                    s.m1   = m1;
                    s.m2   = m2;
            end
            s.m1 = m1;
            s.m2 = m2;
            sM = SmoothingExponentComputer.create(s);
            obj.superEllipseExponent = sM;
        end
        
        function createRotator(obj)
            obj.rotator = ConstitutiveTensorRotator();
        end
        
        function obtainReferenceConsistutiveTensor(obj,x)
            mx = x{1};
            my = x{2};
            [c,dc] = obj.Ctensor.compute([mx,my]);
            obj.Cref  = c;
            obj.dCref = permute(dc,[1 2 4 3 5]);
        end
        
        function obtainReferenceAmplificatorTensor(obj,x,pNorm)
            obj.xV = x;
            mx = x{1};
            my = x{2};
            [P2,dP2,Pp,dPp] = obj.Ptensor.compute([mx,my],pNorm);
            obj.P2Ref  = P2;
            obj.dP2Ref = permute(dP2,[1 2 4 3]);
            obj.PpRef  = Pp;
            obj.dPpRef = dPp;
        end
      
        function rotateConstitutiveTensor(obj,alpha)
            obj.rotator.evaluateRotationMatrixForC(alpha);
            cr  = obj.Cref;
            dCr = obj.dCref;
            
           c           = obj.rotator.rotateC(cr);
           dc(:,:,1,:) = obj.rotator.rotateC(squeeze(dCr(:,:,1,:)));
           dc(:,:,2,:) = obj.rotator.rotateC(squeeze(dCr(:,:,2,:)));

            %c = cr;
            %dc = dCr;
            
            obj.C = c;
            obj.dC = dc;
        end
  
        function rotateAmplificatorTensor(obj,alpha)
            obj.rotator.evaluateRotationMatrixForP(alpha);
            
            
            %obj.Pp           = obj.rotator.rotateP(obj.P2Ref);
            %obj.dPp(:,:,1,:) = obj.rotator.rotateP(squeeze(obj.dP2Ref(:,:,1,:)));
            %obj.dPp(:,:,2,:) = obj.rotator.rotateP(squeeze(obj.dP2Ref(:,:,2,:)));
            
            obj.Pp  = obj.PpRef;
            obj.dPp = obj.dPpRef;
            
        %             obj.Pp(1,:) = 0;
        %             obj.Pp(2,:) = 0;
        %             obj.Pp(3,:) = 0;
        %             obj.Pp(4,:) = 2;
        %             obj.Pp(5,:) = 1;
        %             obj.Pp(6,:) = 1;
        %             obj.dPp(:,:,:) = 0;
        end
        
        function q = computeSuperEllipseExponent(obj,x)
            xR = obj.reshapeDesignVariable(x);
            m1 = xR{1};
            m2 = xR{2};
            q = obj.computeSuperEllipse(m1,m2);
        end
        
        function q = computeSuperEllipse(obj,m1,m2)
            obj.createSuperEllipseExponent(m1,m2);
            q(:,1) = obj.superEllipseExponent.compute();
        end
        
        function computeXiSymmetry(obj,x)
            sRelator = SuperEllipseParamsRelator;
            xR = obj.reshapeDesignVariable(x);
            m1 = xR{1};
            m2 = xR{2};
            [rho,~] = obj.density.compute([m1,m2]);
            xi   = sRelator.xi(m1,m2);
            xiUB = sRelator.xiUB(rho);
            xiS  = abs(xi-pi/4)./(xiUB - pi/4);
            obj.xiSymmetry = 1 - xiS;
        end
        
    end

    
end