classdef ShFunc_ComplianceComparison_constraint < ShFunWithElasticPdes
    
    properties (Access = private)
        compliance
        fieldToPrint
        adjointProblem
        comparisonComplianceFrac = 4.7
    end
    
    methods (Access = public)
        
        function obj = ShFunc_ComplianceComparison_constraint(cParams)
            cParams.filterParams.quadratureOrder = 'LINEAR';
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
                        'VectorGauss'};
            names = {'Primal','ComplianceGauss','AlphaGauss',...
                        'AlphaAbsGauss'};
            fP = obj.obtainPrintVariables(types,names);
            fP = obj.addHomogPrintVariablesNames(fP);
        end
        
    end
    
    methods (Access = protected)
        
        function solveState(obj)
            obj.physicalProblem.setC(obj.homogenizedVariablesComputer.C) % (:,:,7200,4); cmat
%             obj.physicalProblem.computeStiffnessMatrix();
%             obj.physicalProblem.computeVariables();
            obj.physicalProblem.solve();
        end
        
        function solveAdjoint(obj)
            obj.adjointProblem = obj.physicalProblem;
        end
        
        function computeFunctionValue(obj)
            phy = obj.physicalProblem;
            dvolum = phy.getDvolume()';
            stress = phy.variables.stress;
            strain = phy.variables.strain;
            ngaus  = size(strain,1);
            nelem  = size(strain,3);

            c = zeros(nelem,ngaus);
            for igaus = 1:ngaus
                stressG = squeeze(stress(igaus,:,:));
                strainG  = squeeze(strain(igaus,:,:));
                e = stressG.*strainG;
                c(:,igaus) = c(:,igaus) + sum(e)';
            end
            obj.compliance = c;
            int = c.*dvolum;
            obj.value = sum(int(:));
            if isempty(obj.value0)
                obj.value = obj.value;
            else
                obj.value = obj.value - obj.comparisonComplianceFrac*obj.value0;
            end
        end
        
        function computeGradientValue(obj)
            phy = obj.physicalProblem;
            ep    = phy.variables.strain;
            ngaus  = size(ep,1);
            nstre  = size(ep,2);
            nelem  = size(ep,3);
            g = zeros(nelem,ngaus,obj.nVariables);
            for igaus = 1:ngaus
                for istre = 1:nstre
                    for jstre = 1:nstre
                        eu_i = squeeze(ep(igaus,istre,:));
                        ep_j = squeeze(ep(igaus,jstre,:));
                        for ivar = 1:obj.nVariables
                            dCij = squeeze(obj.homogenizedVariablesComputer.dC(istre,jstre,ivar,:,igaus));
                            g(:,igaus,ivar) = g(:,igaus,ivar) + (-eu_i.*dCij.*ep_j);
                        end
                    end
                end
            end
            obj.gradient = g;
        end
        
        function f = getPdesVariablesToPrint(obj)
            f{1} = obj.getPdeVariableToPrint(obj.physicalProblem);
        end
        
    end
    
end

