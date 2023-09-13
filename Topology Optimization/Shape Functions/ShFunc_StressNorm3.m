classdef ShFunc_StressNorm3 < ShFunWithElasticPdes
    
    properties (Access = private)
        adjointProb
        fAdjoint
        stressNorm
        stressNormSquare
        stressNorml2LpNorm
        isElementToOptimize
        integralSigmaP
        
        ngaus
        nelem
        sigma
        sigmaA
        sigma2
        sigmaP
        
        
        dsigHPdsigH2v
        PpdSigHP_DsigH2M
    end
    
    methods (Access = public)
        
        function obj = ShFunc_StressNorm3(cParams)
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
        end
        
        function f = getPdesVariablesToPrint(obj)
            f{1} = obj.getPdeVariableToPrint(obj.physicalProblem);
            f{2} = obj.getPdeVariableToPrint(obj.adjointProb);
        end
    
       function fP = addPrintableVariables(obj)
            fP{1}.value = obj.getPdeVariableToPrint(obj.physicalProblem);
            fP{2}.value = obj.getPdeVariableToPrint(obj.adjointProb);
            fP{3}.value = sqrt(obj.sigma2');
            fP{4}.value = obj.designVariable.alpha;
            fP{5}.value = abs(obj.designVariable.alpha);
            fP = obj.addHomogVariables(fP);
       end
        
        function fP = createPrintVariables(obj)
            types = {'Elasticity','Elasticity','ScalarGauss'...
                        'VectorGauss','VectorGauss'};
            names = {'Primal','Adjoint','StressNormGauss',...
                        'AlphaGauss','AlphaAbsGauss'};
            fP = obj.obtainPrintVariables(types,names);
            fP = obj.addHomogPrintVariablesNames(fP);
        end
        
        function v = getVariablesToPlot(obj)
            v{1} = max(sqrt(obj.sigma2));
            sNorm2 = obj.computeClosedSigmaProduct(obj.sigma,obj.sigma);
            v{2} = max(sqrt(sNorm2));
        end
        
        function t = getTitlesToPlot(obj)
            t{1} = 'Max amplified stress';
            t{2} = 'Max stress';
        end
        
        function setPnorm(obj,pNorm)
            obj.target_parameters.stressNormExponent = pNorm;
        end
        
    end
    
    methods (Access = protected)
        
        function updateHomogenizedMaterialProperties(obj)
            obj.filterDesignVariable();
            obj.homogenizedVariablesComputer.computeCtensor(obj.regDesignVariable);
            pNorm = obj.obtainPnorm();
            obj.homogenizedVariablesComputer.computePtensor(obj.regDesignVariable,pNorm);
        end
        
        function computeFunctionValue(obj)
            p  = obj.obtainPnorm(); 
            Pp = obj.homogenizedVariablesComputer.Pp;
            obj.sigma  = obj.computeSigmaTensor();
            obj.sigmaA = obj.amplifySigma(Pp,obj.sigma);
            obj.sigma2 = obj.computeClosedSigmaProduct(obj.sigmaA,obj.sigmaA);
            obj.sigmaP = obj.computeSigmaP(obj.sigma2);
            obj.integralSigmaP = obj.integrateSigmaP();
            obj.value = obj.integralSigmaP.^(1/p);
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
            obj.computeFadjoint();
            obj.adjointProb.setC(obj.homogenizedVariablesComputer.C);
            obj.adjointProb.computeStiffnessMatrix();
            obj.adjointProb.computeVariablesWithBodyForces(obj.fAdjoint);
        end
        
    end
    
    methods (Access = private)
        
        function g = initG(obj)
            g = zeros(obj.nelem,obj.ngaus,obj.nVariables);
        end
        
        function euR = rotateStrain(obj,eu)
            rot = obj.homogenizedVariablesComputer.rotator;
            euR  = rot.rotateStrain(eu);
        end
        
        function sR = rotateStress(obj,s)
            rot = obj.homogenizedVariablesComputer.rotator; 
            sR  = rot.rotateStress(s);
        end
        
        function ds2ds = computedSigma2dSigma(obj,ds)
            Pp    = obj.homogenizedVariablesComputer.Pp;
            s     = obj.sigma;
            sA    = obj.amplifySigma(Pp,s);
            dsA   = obj.amplifySigma(Pp,ds);
            sXds  = obj.computeClosedSigmaProduct(sA,dsA);
            dsXs  = obj.computeClosedSigmaProduct(dsA,sA);
            ds2ds = sXds + dsXs;
        end
        
        function ds2dP = computedSigma2dP(obj,dPp)
            s     = obj.sigma;
            Pp    = obj.homogenizedVariablesComputer.Pp;
            sA    = obj.amplifySigma(Pp,s);
            dsA   = obj.amplifySigma(dPp,s);
            sXds  = obj.computeClosedSigmaProduct(sA,dsA);
            dsXs  = obj.computeClosedSigmaProduct(dsA,sA);
            ds2dP = sXds + dsXs;
        end
    
        function s = computeSigmaTensor(obj)
            eu = obj.physicalProblem.variables.strain;
            C  = obj.homogenizedVariablesComputer.Cref;
            s  =  obj.computeSigma(C,eu);
        end
        
        function sA = amplifySigma(obj,Pp,s)
            sA = obj.computeSigma(Pp,s);
        end
        
        function computeSigma2(obj)
            s = obj.sigmaA;
            obj.sigma2  = obj.computeClosedSigmaProduct(s,s);
        end
        
        function sNorm = computeClosedSigmaProduct(obj,sa,sb)
            nStre = size(sa,2);
            factor = [1 1 2];
            sNorm = zeros(obj.ngaus,obj.nelem);
            for iStre = 1:nStre
                sAi = obj.squeezeParticular(sa(:,iStre,:),2);
                sBi = obj.squeezeParticular(sb(:,iStre,:),2);
                sNorm = sNorm + factor(iStre)*(sAi.*sBi);
            end
        end
        
        function intSigmaP = integrateSigmaP(obj)
            dvolum = obj.physicalProblem.geometry.dvolu';
            int    = obj.sigmaP.*dvolum;
            intSigmaP = int(:,obj.isElementToOptimize);
            intSigmaP = sum(intSigmaP(:));
        end
        
        function sP = computeSigmaP(obj,s2)
            p  = obj.obtainPnorm;
            sP = s2.^(p/2);
        end
        
        function computeFadjoint(obj)
            dSigmaPdu    = obj.computedSigmaPdu();
            dJdu         = obj.computedJdSigmaP(dSigmaPdu);
            obj.fAdjoint = obj.assambleVector(-dJdu);
        end
        
        function dJdx = computeDJdx(obj)
            dSigmaPdx = obj.computedSigmaPdx();
            dJdx      = obj.computedJdSigmaP(dSigmaPdx);
            dJdx(~obj.isElementToOptimize,:) = 0;
        end

        function s = computeSigma(obj,C,eu)
            %eu = obj.rotateStrain(eu);
            s  = obj.computeStress(C,eu);
            %s = obj.rotateStress(s);
        end
        
        function dSigmaPdx = computedSigmaPdx(obj)
            eu = obj.physicalProblem.variables.strain;
            dC = obj.homogenizedVariablesComputer.dCref;
            dP = obj.homogenizedVariablesComputer.dPp;
            nX = size(dC,3);
            dSigmaPdx = zeros(nX,obj.ngaus,obj.nelem);
            for ix = 1:nX
                dCdx       = squeeze(dC(:,:,ix,:));
                dPdx       = squeeze(dP(:,:,ix,:));
                dSigmadx   = obj.computeSigma(dCdx,eu);
                dSigma2dxA = obj.computedSigma2dSigma(dSigmadx); 
                dSigma2dxB = obj.computedSigma2dP(dPdx); 
                dSigma2dx  = dSigma2dxA + dSigma2dxB;
                dSigmaPdx(ix,:,:) = obj.computedSigmaPdSigma2(dSigma2dx);
            end
            dSigmaPdx = permute(dSigmaPdx,[3 2 1]);
        end
        
        function dSigmaPdu = computedSigmaPdu(obj)
            C    = obj.homogenizedVariablesComputer.Cref;
            dEps = obj.computedEps();
            nV   = size(dEps,1);
            dSigmaPdu = zeros(nV,obj.ngaus,obj.nelem);
            for iv = 1:nV
                dEpsDu    = obj.squeezeParticular(dEps(iv,:,:,:),1);
                dSigmadu  = obj.computeSigma(C,dEpsDu);
                dSigma2du = obj.computedSigma2dSigma(dSigmadu);
                dSigmaPdu(iv,:,:) = obj.computedSigmaPdSigma2(dSigma2du);
            end
        end
        
        function dSigmaPdx = computedSigmaPdSigma2(obj,dSigma2dx)
            p = obj.obtainPnorm();
            dSigmaPdSigma2 = p/2*(obj.sigma2).^(p/2-1);
            dSigmaPdx = dSigmaPdSigma2.*dSigma2dx;
        end

        function s = computeStress(obj,C,e)
            s = zeros(size(e));
            for igaus = 1:obj.ngaus
                eG(1,:,:) = squeeze(e(igaus,:,:));
                sGt = bsxfun(@times,C,eG);
                s(igaus,:,:) = obj.squeezeParticular(sum(sGt,2),2);
            end
        end
        
        function g = computeDHdx(obj)
            eu = obj.physicalProblem.variables.strain;
            %eu = obj.rotateStrain(eu);
            dC = obj.homogenizedVariablesComputer.dCref;
            ep = obj.adjointProb.variables.strain;
            g = obj.initG();
            for ivar = 1:obj.nVariables
                dCiv = squeeze(dC(:,:,ivar,:));
                ds = obj.computeSigma(dCiv,ep);
                g(:,:,ivar) = obj.squeezeParticular(sum(eu.*ds,2),2)';
            end
        end
        
        function p = obtainPnorm(obj)
            p = obj.target_parameters.stressNormExponent;
        end
        
        function createAdjointProblem(obj,fileName)
            obj.adjointProb = FEM.create(fileName);
        end
         
        function dEps = computedEps(obj)
            dvolum = obj.physicalProblem.geometry.dvolu';
            nstre = obj.physicalProblem.element.getNstre();
            nnode = obj.physicalProblem.element.nnode;
            nunkn = obj.physicalProblem.element.dof.nunkn;
            nv    = nnode*nunkn;
            dEps = zeros(nv,obj.ngaus,nstre,obj.nelem);
            for igaus = 1:obj.ngaus
               B  = obj.physicalProblem.element.computeB(igaus);
               Bm = permute(B,[2 1 3]);
               dvG(1,1,:) = squeeze(dvolum(igaus,:));
               dvGm = repmat(dvG,nv,nstre,1);
               dEps(:,igaus,:,:) = Bm.*dvGm;
            end
        end
        
        function fAdjoint = assambleVector(obj,dsigmaPdu)
            eforce = dsigmaPdu;
            eforce(:,:,~obj.isElementToOptimize) = 0;
            Fvol = obj.physicalProblem.element.AssembleVector({eforce});
            fAdjoint = Fvol;
        end
        
        function dJdx = computedJdSigmaP(obj,dSigmaPdx)
            p = obj.obtainPnorm();
            intSigmaP = obj.integralSigmaP;
            dJ_dSigmaP = 1/p*intSigmaP^(1/p-1);
            dJdx = dJ_dSigmaP*dSigmaPdx;
        end
        
        function createElementsToOptimize(obj)
             phy = obj.getPhysicalProblems();
             m = phy{1}.mesh;
             xG = transpose(m.computeBaricenter);
             x = xG(:,1);
             y = xG(:,2);
%             isDirichletPartX = x > -1e-12 & x < 0.05;
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
        
       
    end
    
    methods (Access = private, Static)
    
        function squeezed = squeezeParticular(t,dim)
            z = size(t);
            index = setdiff(1:length(z),dim);
            squeezed = reshape(t,[z(index) z(dim)]);
        end
    
    end
    
end
