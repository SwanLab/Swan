classdef ShFunc_Chomog < ShapeFunctional
    
    properties (Access = public)
        Chomog
        tstress
        tstrain
        dCh
        physicalProblem
        rho
        matProps
    end
    
    properties (Access = protected)
        regDesignVariable
    end
    
    methods (Access = public)
        
        
        function f = getPhysicalProblems(obj)
            f{1} = obj.physicalProblem;
        end
        
        function f = getRegularizedDensity(obj)
            f = obj.regDesignVariable;
        end
        
        function fP = addPrintableVariables(obj)
            fP = obj.addMicroVariables();
            fP{end+1}.value = obj.gradient;
            fH = obj.homogenizedVariablesComputer.addPrintableVariables(obj.designVariable);
            for i = 1:numel(fH)
                fP{end+1} = fH{i};
            end
        end
        
        function fP = addMicroVariables(obj)
            p = obj.getPhysicalProblems();
            p = p{1};
            nMicro = numel(p.variables2print);
            fP = cell(nMicro,1);
            for i = 1:nMicro
                pV = p;
                pV.variables = p.variables2print{i};
                cParams.physicalProblem = pV;
                g = PdeVariableToPrintGetter(cParams);
                v = g.compute();
                fP{i}.value = v;
            end
        end
        
        function fP = createPrintVariables(obj)
            types = {'Elasticity','Elasticity','Elasticity','ScalarNodal'};
            names = {'Primal1','Primal2','Primal3','Gradient'};
            
            fP = obj.obtainPrintVariables(types,names);
            fP = obj.addHomogPrintVariablesNames(fP);
        end
        
        function v = getVariablesToPlot(obj)
            v{1} = obj.value*obj.value0;
        end
        
        function t = getTitlesToPlot(obj)
            t{1} = 'C - C not scaled';
        end
        
        function q = getQuad(obj)
            q = obj.physicalProblem.element.quadrature;
        end
        
        function computeFunction(obj)
            obj.updateHomogenizedMaterialProperties();
            obj.solveState();
            obj.computeFunctionValue();
            obj.normalizeFunction();
        end
        
        function computeGradient(obj)
            obj.computeGradientValue();
            obj.filterGradient();
            obj.normalizeGradient();
        end
        
        function computeFunctionAndGradient(obj)
            obj.computeFunction();
            obj.computeGradient();
        end
        
    end
    
    methods (Access = protected)
        
        function filterGradient(obj)
            g = obj.gradient;
            gf = zeros(size(obj.Msmooth,1),obj.nVariables);
            for ivar = 1:obj.nVariables
                gs = g(:,:,ivar);
                gf(:,ivar) = obj.filter.getP1fromP0(gs);
            end
            %gf = obj.Msmooth*gf;
            g = gf(:);
            obj.gradient = g;
        end                
        
        function computeChDerivative(obj)
            dC = obj.homogenizedVariablesComputer.dC;
            
            nStre = obj.physicalProblem.element.getNstre();
            nelem = obj.physicalProblem.mesh.nelem;
            ngaus = obj.physicalProblem.element.quadrature.ngaus;
            
            
            dChV = zeros(nStre,nStre,nelem,ngaus);
            for iStre = 1:nStre
                for jStre = 1:nStre
                    for igaus=1:ngaus
                        for kStre =1:nStre
                            for lStre = 1:nStre
                                dChVij = squeeze(dChV(iStre,jStre,:,igaus));
                                eik   = squeeze(obj.tstrain(iStre,igaus,kStre,:));
                                dCkl  = squeeze(dC(kStre,lStre,:,:,igaus));
                                ejl   = squeeze(obj.tstrain(jStre,igaus,lStre,:));
                                dChV(iStre,jStre,:,igaus) = dChVij + eik.*dCkl.*ejl;%Original
                            end
                        end
                    end
                end
            end
            obj.dCh = dChV;
        end
        
        function initChomog(obj,cParams)
            cParams.filterParams.quadratureOrder = 'LINEAR';
            obj.init(cParams);
            obj.physicalProblem = FEM.create(cParams.femSettings.fileName);
        end
        
        function solveState(obj)
            obj.physicalProblem.setC(obj.homogenizedVariablesComputer.C)
            obj.physicalProblem.computeChomog();
            obj.Chomog  = obj.physicalProblem.variables.Chomog;
            obj.tstrain = obj.physicalProblem.variables.tstrain;
            obj.tstress = obj.physicalProblem.variables.tstress;
        end
        
        function updateHomogenizedMaterialProperties(obj)
            rhoV{1} = obj.filter.getP0fromP1(obj.designVariable.value);
            obj.regDesignVariable = rhoV{1};
            obj.homogenizedVariablesComputer.computeCtensor(rhoV);
            obj.homogenizedVariablesComputer.computeDensity(rhoV);
        end
        
        
        function dChInv = computedChInv(obj,Ch,alpha,beta)
            invCh      = inv(Ch);
            nStres     = obj.getnStre();
            weights    = alpha*beta';
            weightsInv = invCh*weights*invCh;
            dChInv     = zeros(size(obj.dCh,[3 4]));
            for iStre = 1:nStres
                for jStre = 1:nStres
                    DtCij  = squeezeParticular(obj.dCh(iStre,jStre,:,:),[1 2]);
                    wij    = weightsInv(iStre,jStre);
                    dChInv = dChInv - wij*DtCij;
                end
            end
        end
        
        function n = getnStre(obj)
            n = obj.physicalProblem.element.getNstre();
        end
        
        function n = getnElem(obj)
            n = obj.physicalProblem.element.nelem;
        end
        
        function n = getnGaus(obj)
            n = size(obj.tstrain,2);
        end
        
    end
    
    methods (Access = protected, Static)
        
        function aCb = projectTensor(C,a,b)
            ab  = a*b';
            aCb = ab.*C;
            aCb = sum(aCb(:));
        end
        
    end
end
