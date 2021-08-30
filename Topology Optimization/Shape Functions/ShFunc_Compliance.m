classdef ShFunc_Compliance < ShFunWithElasticPdes
    
    properties (Access = private)
        compliance
        fieldToPrint
        adjointProblem
    end
    
    methods (Access = public)
        
        function obj = ShFunc_Compliance(cParams)
            cParams.filterParams.quadratureOrder = 'LINEAR';
            obj.init(cParams);
            fileName = cParams.femSettings.fileName;
            obj.createEquilibriumProblem(fileName);
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
            obj.physicalProblem.setC(obj.homogenizedVariablesComputer.C)
            obj.physicalProblem.computeVariables();
        end
        
        function solveAdjoint(obj)
            obj.adjointProblem = obj.physicalProblem;
        end
        
        function computeFunctionValue(obj)
            phy = obj.physicalProblem;
            dvolum  = phy.geometry.dvolu;
            stress = phy.variables.stress;
            strain = phy.variables.strain;            
            ngaus  = phy.element.quadrature.ngaus;
            nelem  = phy.mesh.nelem;
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
        end
                    
        function computeGradientValue(obj)
            eu    = obj.physicalProblem.variables.strain;
            ep    = obj.physicalProblem.variables.strain;
            nelem = obj.physicalProblem.mesh.nelem;
            ngaus = obj.physicalProblem.element.quadrature.ngaus;
            nstre = obj.physicalProblem.element.getNstre();             
            g = zeros(nelem,ngaus,obj.nVariables);
            for igaus = 1:ngaus
                for istre = 1:nstre
                    for jstre = 1:nstre
                        eu_i = squeeze(eu(igaus,istre,:));
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

