classdef ShFunc_Compliance < ShFunWithElasticPdes
    
    properties (Access = private)
        compliance
        fieldToPrint
        adjointProblem
        gradientGauss
    end
    
    methods (Access = public)
        
        function obj = ShFunc_Compliance(cParams)
            obj.init(cParams);
            obj.physicalProblem = cParams.femSettings.physicalProblem;
            obj.createOrientationUpdater();
        end
        
        function fP = addPrintableVariables(obj)
            phy = obj.getPdesVariablesToPrint();
            fP{1}.value = phy{1};
            fP{2}.value = obj.compliance/obj.value0;
            fP{3}.value = obj.designVariable.alpha;
            fP{4}.value = abs(obj.designVariable.alpha);
            fP{5}.value = permute(obj.gradientGauss,[3 1 2]);
            fP = obj.addHomogVariables(fP);
        end

        function v = getVariablesToPlot(obj)
            v{1} = obj.value*obj.value0;
        end
        
        function t = getTitlesToPlot(obj)
            t{1} = 'Compliance non scaled';
        end
        
        function fP = createPrintVariables(obj)
            types = {'Elasticity','ScalarGauss','VectorGauss'...
                        'VectorGauss','VectorGauss'};
            names = {'Primal','ComplianceGauss','AlphaGauss',...
                        'AlphaAbsGauss','GradientGauss'};
            fP = obj.obtainPrintVariables(types,names);
            fP = obj.addHomogPrintVariablesNames(fP);
        end
        
        function [fun, funNames] = getFunsToPlot(obj)
            mesh = obj.designVariable.mesh;
            phy = obj.physicalProblem;
            strain = phy.strainFun; % !!!
            stress = phy.stressFun; % !!!
            displ  = phy.uFun; % !!!
            compl  = obj.compliance/obj.value0;

            quad = Quadrature.set(mesh.type);
            quad.computeQuadrature('LINEAR');

            aa.mesh       = mesh;
            aa.quadrature = quad;
            aa.fValues    = permute(compl, [3 2 1]);
            complFun = FGaussDiscontinuousFunction(aa);
            
            bb.mesh    = mesh;
            bb.fValues = obj.designVariable.alpha';
            alphaFun = P0Function(bb);

            fun      = {complFun, strain, stress, displ};
            funNames = {'compliance', 'strain', 'stress', 'u'};
        end
    end
    
    methods (Access = protected)
        
        function solveState(obj)
            obj.physicalProblem.setC(obj.homogenizedVariablesComputer.C) % (:,:,7200,4); cmat
            obj.physicalProblem.solve();
        end
        
        function solveAdjoint(obj)
            obj.adjointProblem = obj.physicalProblem;
        end
        
        function computeFunctionValue(obj)
            phy = obj.physicalProblem;
            stress = phy.stressFun.fValues;
            strain = phy.strainFun.fValues;
            dvolum = obj.designVariable.mesh.computeDvolume(phy.strainFun.quadrature)';
            ngaus  = size(strain,2);
            nelem  = size(strain,3);

            c = zeros(nelem,ngaus);
            for igaus = 1:ngaus
                stressG = squeeze(stress(:,igaus,:));
                strainG  = squeeze(strain(:,igaus,:));
                e = stressG.*strainG;
                c(:,igaus) = c(:,igaus) + sum(e)';
            end
            obj.compliance = c;
            int = c.*dvolum;
            obj.value = sum(int(:));
        end
        
        function computeGradientValue(obj)
            obj.computeGradientInGauss();
            obj.gradient = obj.gradientGauss;
        end

        function computeGradientInGauss(obj)
            phy = obj.physicalProblem;
            ep    = phy.strainFun.fValues;
            nstre  = size(ep,1);
            ngaus  = size(ep,2);
            nelem  = size(ep,3);
            g = zeros(nelem,ngaus,obj.nVariables);
            for igaus = 1:ngaus
                for istre = 1:nstre
                    for jstre = 1:nstre
                        eu_i = squeeze(ep(istre,igaus,:));
                        ep_j = squeeze(ep(jstre,igaus,:));
                        for ivar = 1:obj.nVariables
                            dCij = squeeze(obj.homogenizedVariablesComputer.dC(istre,jstre,ivar,:,igaus));
                            g(:,igaus,ivar) = g(:,igaus,ivar) + (-eu_i.*dCij.*ep_j);
                        end
                    end
                end
            end            
            obj.gradientGauss = g;
        end
        
        function f = getPdesVariablesToPrint(obj)
            f{1} = obj.getPdeVariableToPrint(obj.physicalProblem);
        end
        
    end
    
end

