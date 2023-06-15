classdef ShFunc_StressNorm2 < ShFunWithElasticPdes
    
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
        sigma2
        sigmaP
        
        
        dsigHPdsigH2v
        PpdSigHP_DsigH2M
    end
    
    methods (Access = public)
        
        function obj = ShFunc_StressNorm2(cParams)
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
        
        function setPnorm(obj,pNorm)
            obj.target_parameters.stressNormExponent = pNorm;            
        end
        
    end
    
    methods (Access = protected)
        
        function updateHomogenizedMaterialProperties(obj)
            obj.filterDesignVariable();
            obj.homogenizedVariablesComputer.computeCtensor(obj.regDesignVariable);
            pNormIndex = obj.computePNormIndex();
            obj.homogenizedVariablesComputer.computePtensor(obj.regDesignVariable,pNormIndex);
        end
        
        function computeFunctionValue(obj)          
            p = obj.obtainPnorm(); 
            obj.computeSigmaTensor();
            obj.computeSigma2();
            obj.computeSigmaP();
            obj.integrateSigmaP();            
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
        
        function ds2M = computedSigma2dSigma(obj,s,ds)
            sXds = obj.computeClosedSigmaProduct(s,ds);
            dsXs = obj.computeClosedSigmaProduct(ds,s);
            ds2M = sXds + dsXs;            
        end
        
        function computeSigmaTensor(obj)
            eu = obj.physicalProblem.variables.strain;
            C  = obj.homogenizedVariablesComputer.Cref;   
            obj.sigma  = obj.computeSigma(C,eu);        
        end
        
        function computeSigma2(obj)
            s = obj.sigma;
            obj.sigma2  = obj.computeClosedSigmaProduct(s,s);
        end
        
        function sNorm = computeClosedSigmaProduct(obj,sa,sb)
            nStre = size(sa,2);
            factor = [1 1 2];
            sNorm = zeros(obj.ngaus,1,obj.nelem);
            for iStre = 1:nStre
                sAi = sa(:,iStre,:);                
                sBi = sb(:,iStre,:);
                sNorm = sNorm + factor(iStre)*(sAi.*sBi);
            end            
            sNorm = obj.squeezeParticular(sNorm,1);
        end
        
        function integrateSigmaP(obj)          
            dvolum = obj.physicalProblem.geometry.dvolu';          
            int    = obj.sigmaP.*dvolum;            
            intOpt = int(:,obj.isElementToOptimize);
            obj.integralSigmaP = sum(intOpt(:));
        end        
        
        function computeSigmaP(obj)
            s2 = obj.sigma2;
            p  = obj.obtainPnorm;
            obj.sigmaP = s2.^(p/2);                          
        end
        
        function computeFadjoint(obj)   
            dSigmaPdu  = obj.computedSigmaPdu();
            dJ_dSigmaP = obj.computedJdSigmaP();
            dSigmaPdu  = -dJ_dSigmaP*dSigmaPdu;            
            obj.fAdjoint = obj.assambleVector(dSigmaPdu);           
        end         
        
        function g = computeDJdx(obj)
            dSigmaPdx  = obj.computedSigmaPdx();
            dJ_dSigmaP = obj.computedJdSigmaP();            
            g = dJ_dSigmaP*dSigmaPdx;            
            g(~obj.isElementToOptimize,:) = 0;
        end

        function s = computeSigma(obj,C,eu)
            eu = obj.rotateStrain(eu);            
            s  = obj.computeStress(C,eu);
            %s = obj.rotateStress(s);                        
        end
        
        function dSigmaPdx = computedSigmaPdx(obj)
            eu = obj.physicalProblem.variables.strain;
            dC = obj.homogenizedVariablesComputer.dCref;
            nX = size(dC,3);       
            dSigmaPdx = zeros(nX,obj.ngaus,obj.nelem);
            for ix = 1:nX
                dCdx      = squeeze(dC(:,:,ix,:));
                dSigmadx  = obj.computeSigma(dCdx,eu);
                dSigma2dx = obj.computedSigma2dSigma(obj.sigma,dSigmadx);                
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
                dSigma2du = obj.computedSigma2dSigma(obj.sigma,dSigmadu);                
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
            eu = obj.rotateStrain(eu);
            dC = obj.homogenizedVariablesComputer.dCref;
            ep = obj.adjointProb.variables.strain;                        
            g = obj.initG();            
            for ivar = 1:obj.nVariables                
                dCiv = squeeze(dC(:,:,ivar,:));
                ds = obj.computeSigma(dCiv,ep);                  
                g(:,:,ivar) = obj.squeezeParticular(sum(eu.*ds,2),2);
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
               dvG = repmat(dvG,nv,nstre,1);
               dEps(:,igaus,:,:) = Bm.*dvG;
            end
        end
        
        function fAdjoint = assambleVector(obj,dsigmaPdu)
            eforce = dsigmaPdu;
            eforce(:,:,~obj.isElementToOptimize) = 0;
            Fvol = obj.physicalProblem.element.AssembleVector({eforce});
            fAdjoint = Fvol;
        end
        
        function v = computedJdSigmaP(obj)
            p = obj.obtainPnorm();            
            intSigmaP = obj.integralSigmaP;
            v = 1/p*intSigmaP^(1/p-1);            
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
        
        function pNormIndex = computePNormIndex(obj)
           p = obj.obtainPnorm();
           pC = 2:2:32;
           isIndex = p == pC;
           pNormIndex = find(isIndex);
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
