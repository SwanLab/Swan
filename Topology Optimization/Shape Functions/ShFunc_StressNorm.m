classdef ShFunc_StressNorm < ShFunWithElasticPdes
    
    properties (Access = private)
        adjointProb
        fAdjoint
    
        isElementToOptimize
        integralSigmaAP
        sigmaAP
        sigmaAPMax
        pMax
        ngaus
        nelem
        
        sigma
        PpdSigP_Dsig2M
       
    end
    
    methods (Access = public)
        
        function obj = ShFunc_StressNorm(cParams)
            cParams.filterParams.quadratureOrder = 'LINEAR';
            obj.init(cParams);
            obj.physicalProblem = cParams.femSettings.physicalProblem;
            fileName = cParams.femSettings.fileName;
            obj.createAdjointProblem(fileName);
            obj.createOrientationUpdater();
            obj.createElementsToOptimize();
            obj.updateHomogenizedMaterialProperties();
            obj.nelem = obj.physicalProblem.mesh.nelem;
            obj.ngaus = obj.physicalProblem.element.quadrature.ngaus;
            obj.pMax = 16;
        end
        
        function f = getPdesVariablesToPrint(obj)
            f{1} = obj.getPdeVariableToPrint(obj.physicalProblem);
            f{2} = obj.getPdeVariableToPrint(obj.adjointProb);
        end
    
       function fP = addPrintableVariables(obj)
            fP{1}.value = obj.getPdeVariableToPrint(obj.physicalProblem);
            fP{2}.value = obj.getPdeVariableToPrint(obj.adjointProb);
            fP{3}.value = (obj.sigmaAP').^(1/obj.obtainPnorm());
            fP{4}.value = obj.designVariable.alpha;
            fP{5}.value = abs(obj.designVariable.alpha);
           % obj.computeSigmaAPMax();
            fP{6}.value = (obj.sigmaAPMax').^(1/obj.pMax);
            fP = obj.addHomogVariables(fP);
       end
        
        function fP = createPrintVariables(obj)
            types = {'Elasticity','Elasticity','ScalarGauss'...
                        'VectorGauss','VectorGauss','ScalarGauss'};
            names = {'Primal','Adjoint','AmplifiedStressNormP',...
                        'AlphaGauss','AlphaAbsGauss','AmplifiedStressNormMax'};
            fP = obj.obtainPrintVariables(types,names);
            fP = obj.addHomogPrintVariablesNames(fP);
        end
        
        function setPnorm(obj,pNorm)
            obj.target_parameters.stressNormExponent = pNorm;
        end

        function v = getVariablesToPlot(obj)
          obj.computeSigmaAPMax();
            intSigmaAMax = obj.integrateSigmaAP(obj.sigmaAPMax);
            sigmaAMax = intSigmaAMax.^(1/obj.pMax);
       %sigmaAMax = 1;
            v{1} = sigmaAMax;%/obj.value0;
            v{2} = sigmaAMax/obj.value0;
            v{3} = obj.computeMaxSigmaNorm();
            v{4} = obj.obtainPnorm;
        end
        
        function t = getTitlesToPlot(obj)
            t{1} = 'Max amplified stress (real)';
            t{2} = 'Max amplified stress (adim)';
            t{3} = 'Max stress';
            t{4} = 'p norm';
        end
        
    end
    
    methods (Access = protected)
        
        function sNorm = computeSigmaNorm(obj)
            sNorm2 = obj.computeClosedSigmaProduct(obj.sigma,obj.sigma);
            sNorm = sqrt(sNorm2);
        end
        
        function sNormMax = computeMaxSigmaNorm(obj)
            sNorm = obj.computeSigmaNorm();
            sNormMax = max(sNorm);
        end
        
        function computeSigmaAPMax(obj)
            p    = obj.obtainPnorm();
            obj.computeFunctionValueForCertainP(obj.pMax);
            obj.sigmaAPMax = obj.sigmaAP;
            obj.computeFunctionValueForCertainP(p);
        end
        
        function computeFunctionValueForCertainP(obj,p)
            obj.target_parameters.stressNormExponent = p;
            obj.updateHomogenizedMaterialProperties();
            obj.computeFunctionValue();
            obj.normalizeFunction();
            obj.updateAlpha();
        end
        
        function updateHomogenizedMaterialProperties(obj)
            obj.filterDesignVariable();
            obj.homogenizedVariablesComputer.computeCtensor(obj.regDesignVariable);
            pNorm = obj.obtainPnorm();
            obj.homogenizedVariablesComputer.computePtensor(obj.regDesignVariable,pNorm);
        end
        
        function computeFunctionValue(obj)
            p = obj.obtainPnorm();
            obj.sigma   = obj.computeSigmaTensor();
            obj.sigmaAP = obj.computeSigmaAP();
            obj.integralSigmaAP = obj.integrateSigmaAP(obj.sigmaAP);
            obj.value = obj.integralSigmaAP.^(1/p);
        end
        
        function s = computeSigmaTensor(obj)
            eu = obj.physicalProblem.variables.strain;
            C  = obj.homogenizedVariablesComputer.Cref;
            s  = obj.computeSigma(C,eu);
        end
        function solveState(obj)
            obj.physicalProblem.setC(obj.homogenizedVariablesComputer.C);
            obj.physicalProblem.computeStiffnessMatrix();
            obj.physicalProblem.computeVariables();
        end
        
        function computeGradientValue(obj)
            dJdx = obj.computeDJdx();
            dHdx = obj.computeDHdx();
            obj.gradient = dJdx + dHdx;
        end
        
        function solveAdjoint(obj)
            obj.PpdSigP_Dsig2M = obj.computePpdSigP_Dsig2M(obj.sigma);
            obj.computeFadjoint();
            obj.adjointProb.setC(obj.homogenizedVariablesComputer.C);
            obj.adjointProblem.computeStiffnessMatrix();
            obj.adjointProb.computeVariablesWithBodyForces(obj.fAdjoint);
        end
        
    end
    
    methods (Access = private)
        
        function g = initG(obj)
            g = zeros(obj.nelem,obj.ngaus,obj.nVariables);
        end
        
      function sR = rotateStress(obj,s)
            rot = obj.homogenizedVariablesComputer.rotator; 
            sR  = rot.rotateStress(s);
        end
        
        function euR = rotateStrain(obj,eu)
            rot = obj.homogenizedVariablesComputer.rotator; 
            euR  = rot.rotateStrain(eu);
        end
        
        function ds2M = computeDsigma2M(obj,s,ds)
            sXds = obj.computeOpenProduct(s,ds);
            dsXs = obj.computeOpenProduct(ds,s);
            ds2M = sXds + dsXs;
        end
        
        function s2v = tensorToVector(obj,s2)
            nstre = size(s2,1);
            nQuadStre = nstre*(nstre+1)/2;
            s2v = zeros(nQuadStre,obj.ngaus,obj.nelem);
            factor = [1 1 1 0.5 0.5 0.5];
            for i = 1:nstre
                for j = 1:nstre
                   k = obj.indexTensor2Voigt(i,j);
                   f = factor(k);
                   sij(1,1:obj.ngaus,:) = squeeze(s2(i,j,:,:));
                   s2v(k,:,:) = s2v(k,:,:)+ f*sij;
                end
            end
        end
        
        function s2M = vectorToTensor(obj,s2)
            nstre = 3;
            nt = size(s2,1);
            s2M = zeros(nt,nstre,nstre,obj.ngaus,obj.nelem);
            factor = [1 1 1 0.5 0.5 0.5];
            for i = 1:nstre
                for j = 1:nstre
                  k = obj.indexTensor2Voigt(i,j);
                  f = factor(k);
                  sk = squeezeParticular(s2(:,k,:,:),2);
                  s2M(:,i,j,:,:) = f*sk;
                end
            end
        end
       
        function sP = computeSigmaP(obj,sigH2V,alpha)
            nMonom = size(alpha,1);
            sP = ones(nMonom,obj.ngaus,obj.nelem);
            isZero    = any(alpha<0,2);
            isNotZero = ~isZero;
            sP(isZero,:,:) = 0;
            alphat = alpha(isNotZero,:);
            sigH2Ve(1,:,:,:) = sigH2V;
            t = bsxfun(@power,sigH2Ve,alphat);
            sP(isNotZero,:,:) = squeezeParticular(prod(t,2),2);
        end
        
        function sigmaAP = computeSigmaAP(obj)
            Pp = obj.homogenizedVariablesComputer.Pp;
            alpha = obj.obtainAlpha();
            sigma2M = obj.computeOpenProduct(obj.sigma,obj.sigma);
            sigma2V = obj.tensorToVector(sigma2M);
            sigmaP  = obj.computeSigmaP(sigma2V,alpha);
            sigmaAP = obj.computePpSigmaP(Pp,sigmaP);
        end
        
        function sP = computePpSigmaP(obj,Pp,sigmaP)
            Ppe(:,1,:) = Pp;
            Ppe  = repmat(Ppe,1,obj.ngaus,1);
            sigP = Ppe.*sigmaP;
            sP = squeezeParticular(sum(sigP,1),1);
        end

        function intSAP = integrateSigmaAP(obj,sP)
            dvolum = obj.physicalProblem.geometry.dvolu';
            int    = sP.*dvolum;
            intOpt = int(:,obj.isElementToOptimize);
            intSAP = sum(intOpt(:));
        end
        
      function dJdm = computeDJdx(obj)
            Pp_dSigmaAPdx = obj.computePpdSigmaAPdx();
            dPpdx_SigmaAP = obj.computedPpdxSigmaAP();
            dSigmaAP_dx   = Pp_dSigmaAPdx + dPpdx_SigmaAP;
            dJdm = obj.computedJdSigmaAP(dSigmaAP_dx);
            dJdm(~obj.isElementToOptimize,:) = 0;
      end
        
        function dJdx = computedJdSigmaAP(obj,dSigmaAPdx)
            p = obj.obtainPnorm();
            intSigmaAP = obj.integralSigmaAP;
            dJ_dSigmaP = 1/p*intSigmaAP^(1/p-1);
            dJdx = dJ_dSigmaP*dSigmaAPdx;
        end
      
        function sH = computeSigma(obj,C,eu)
            eu = obj.rotateStrain(eu);
            s  = obj.computeStress(C,eu);
            %sH = obj.rotateStress(s);
            sH = s;
        end
        
        function g = computePpdSigmaAPdx(obj)
            eu = obj.physicalProblem.variables.strain;
            dC = obj.homogenizedVariablesComputer.dCref;
            nM    = size(dC,3);
            g = zeros(nM,obj.ngaus,obj.nelem);
            for im = 1:nM
                dCdx        = squeeze(dC(:,:,im,:));
                dSigmadx    = obj.computeSigma(dCdx,eu);
                PpdSigmaPdx = obj.computePpDSigmaPdx(dSigmadx);
                g(im,:,:)   = PpdSigmaPdx;
            end                
            g = permute(g,[3 2 1]);
        end
        
        function PpdSigmaPdx = computePpDSigmaPdx(obj,dSigmadx)
            dSig2Mdx     = obj.computeDsigma2M(obj.sigma,dSigmadx);
            PpdSigmaPdxE = obj.PpdSigP_Dsig2M.*dSig2Mdx;
            PpdSigmaPdx  = squeezeParticular(sum(PpdSigmaPdxE,[1 2]),[1 2]);
        end
        
        function dSigmaAPdu = computedSigmaAPdu(obj)
            C    = obj.homogenizedVariablesComputer.Cref;
            dEps = obj.physicalProblem.element.dEps;
            nV   = size(dEps,1);
            dSigmaAPdu = zeros(nV,obj.ngaus,obj.nelem);
            for iv = 1:nV
                dEpsDu        = squeezeParticular(dEps(iv,:,:,:),1); 
                dSigmadu      = obj.computeSigma(C,dEpsDu);
                PpdSigmaPdu   = obj.computePpDSigmaPdx(dSigmadu);
                dSigmaAPdu(iv,:,:) = PpdSigmaPdu;
            end
        end
 
        function computeFadjoint(obj)
            dSigmaAPdu   = obj.computedSigmaAPdu();
            dJdu         = obj.computedJdSigmaAP(dSigmaAPdu);
            obj.fAdjoint = obj.assambleVector(-dJdu);
        end
 
        function g = computedPpdxSigmaAP(obj)
            alpha = obj.obtainAlpha();
            sigma2M = obj.computeOpenProduct(obj.sigma,obj.sigma);
            sigma2V = obj.tensorToVector(sigma2M);
            sigmaP  = obj.computeSigmaP(sigma2V,alpha);
            dPp = obj.homogenizedVariablesComputer.dPp;
            nM    = size(dPp,3);
            g = zeros(nM,obj.ngaus,obj.nelem);
            for im = 1:nM
                dPpdx       = dPp(:,:,im);
                dPpdxSigmaP = obj.computePpSigmaP(dPpdx,sigmaP);
                g(im,:,:)   = dPpdxSigmaP;
            end
            g = permute(g,[3 2 1]);
        end
        
        function s = computeStress(obj,C,e)
            s = zeros(size(e));
            for igaus = 1:obj.ngaus
                eG(1,:,:) = squeeze(e(igaus,:,:));
                sGt = bsxfun(@times,C,eG);
                s(igaus,:,:) = squeezeParticular(sum(sGt,2),2);
            end
        end
        
        function g = computeDHdx(obj)
            eu = obj.physicalProblem.variables.strain;
            eu = obj.rotateStrain(eu);
            dC = obj.homogenizedVariablesComputer.dCref;
            ep = obj.adjointProb.variables.strain;
            g = obj.initG();
            for ivar = 1:obj.nVariables
                dCiv = squeeze(dC(:,:,ivar,:));
                ds = obj.computeSigma(dCiv,ep);
                g(:,:,ivar) = squeezeParticular(sum(eu.*ds,2),2);
            end
        end
        
        function p = obtainPnorm(obj)
            p = obj.target_parameters.stressNormExponent;
        end
        
        function createAdjointProblem(obj,fileName)
            obj.adjointProb = FEM.create(fileName);
        end
        
        function fAdjoint = assambleVector(obj,dsigmaPdu)
            eforce = dsigmaPdu;
            eforce(:,:,~obj.isElementToOptimize) = 0;
            Fvol = obj.physicalProblem.element.AssembleVector({eforce});
            fAdjoint = Fvol;
        end
            
        function dsigPdsig2v = computeDsigmaP_Dsigma2V(obj,sigma)
            s2v = obj.computeSigma2V(sigma);
            alpha = obj.obtainAlpha();  
            dsigPdsig2v = zeros(size(alpha,1),size(alpha,2),obj.ngaus,obj.nelem);
            alphaP = zeros(size(alpha));
            for iL = 1:size(alpha,2)
                for k = 1:size(alpha,2)
                    if k == iL
                       delta_km = 1;
                    else 
                        delta_km = 0;
                    end
                    alphaP(:,k) = alpha(:,k)-delta_km;
                end
                sP  = obj.computeSigmaP(s2v,alphaP);
                dsigPdsig2v(:,iL,:,:) = bsxfun(@times,alpha(:,iL),sP);
            end
        end
        
        function s2V = computeSigma2V(obj,s)
            s2M = obj.computeOpenProduct(s,s);
            s2V = obj.tensorToVector(s2M);
        end
         
        function alpha = obtainAlpha(obj)
            nstre = obj.physicalProblem.element.getNstre();
            nQuadStre = nstre*(nstre+1)/2;
            p = obj.obtainPnorm();
            alpha  = multinomial_expand(p/2,nQuadStre);
        end
        
        function PpdSigP_Dsig2M = computePpdSigP_Dsig2M(obj,s)
            Pp  = obj.homogenizedVariablesComputer.Pp;
            dsigPdsig2v  = obj.computeDsigmaP_Dsigma2V(s);
            dSigP_Dsig2M = obj.vectorToTensor(dsigPdsig2v);
            Ppe = repmat(Pp,1,1,1,3,3);
            Ppep = permute(Ppe,[1 4 5 3 2]);
            PpdSigP_Dsig2Me = (Ppep.*dSigP_Dsig2M);
            PpdSigP_Dsig2M  = squeezeParticular(sum(PpdSigP_Dsig2Me,1),1);
        end
    
        function createElementsToOptimize(obj)
             phy = obj.getPhysicalProblems();
             m = phy{1}.mesh;
             xG = transpose(m.computeBaricenter);
             x = xG(:,1);
             y = xG(:,2);
%             isDirichletPartX = x > -1e-12 & x < 0.1;
%             isDirichletPartY = y > 0.20 & y < 0.80;
%             isDirichletPart = isDirichletPartX & isDirichletPartY;
%             isNeumannPartX = x > (1-0.05) & x < (1+1e-12);
%             isNeumannPartY = y > 0.40 & y < 0.60;
%             isNeumannPart = isNeumannPartX & isNeumannPartY;
%             isForOptimizing = ~isDirichletPart & ~isNeumannPart;
%             obj.isElementToOptimize = isForOptimizing;

%             s.connec = m.connec(isForOptimizing,:);
%             s.coord  = m.coord;
%             m2 = Mesh(s);
%             m2.plot();
           obj.isElementToOptimize = true(size(x));
        end
        
        function sNorm = computeClosedSigmaProduct(obj,sa,sb)
            nStre = size(sa,2);
            factor = [1 1 2];
            sNorm = zeros(obj.ngaus,obj.nelem);
            for iStre = 1:nStre
                sAi = squeezeParticular(sa(:,iStre,:),2);
                sBi = squeezeParticular(sb(:,iStre,:),2);
                sNorm = sNorm + factor(iStre)*(sAi.*sBi);
            end
        end
 
    end
    
    methods (Access = private, Static)
    
        function k = indexTensor2Voigt(i,j)
            T = [1 6 5;6 2 4;5 4 3];
            k = T(i,j);
        end
         
        function s1s2 = computeOpenProduct(s1,s2)
            ngaus = size(s1,1);
            nstre = size(s1,2);
            nelem = size(s1,3);
            s1s2 = zeros(nstre,nstre,ngaus,nelem);
            for igaus = 1:ngaus
                s1r(1,:,:) = s1(igaus,:,:);
                s2r(1,:,:) = s2(igaus,:,:);
                s1r = permute(s1r,[2 1 3]);
                s1s2(:,:,igaus,:) = bsxfun(@times,s1r,s2r);
            end
        end
    
    end
    
end
