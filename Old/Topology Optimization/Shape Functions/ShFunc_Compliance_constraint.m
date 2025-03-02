classdef ShFunc_Compliance_constraint < ShFunWithElasticPdes
    
    properties (Access = private)
        compliance
        fieldToPrint
        adjointProblem
        targetConstraint
        loadCase
        name
        isFirstIter
    end
    
    methods (Access = public)
        
        function obj = ShFunc_Compliance_constraint(cParams)
            cParams.filterParams.quadratureOrder = 'LINEAR';
            obj.name                             = cParams.type;
            obj.loadCase.number                  = cParams.shNumber;
            obj.init(cParams);
            obj.physicalProblem = cParams.femSettings.physicalProblem;
            obj.createOrientationUpdater();
            obj.isFirstIter = true;
        end

        function t = getTitlesToPlot(obj)
            t{1} = 'Compliance non scaled';
        end

        function v = getVariablesToPlot(obj)
            v{1} = obj.value*obj.value0;
        end
        
    end

    methods (Access = private)

        function c = computeCompliance(obj,n,nv)
            obj.physicalProblem.boundaryConditions.changeBoundaryConditions(n,nv);
            obj.physicalProblem.setC(obj.homogenizedVariablesComputer.C);
            obj.physicalProblem.computeStiffnessMatrix();
            obj.physicalProblem.solve();
            c = obj.computeInitialCompliance();
        end

    end
    
    methods (Access = private)

        function computeTargetConstraint(obj)
            switch obj.name
                case 'complianceConstraintC1'
                    obj.computeMaxComplianceCone();
                    obj.saveLoadCaseCone();
                case 'complianceConstraintC2'
                    obj.computeMaxComplianceCtwo();
                    obj.saveLoadCaseCtwo();
                case 'complianceConstraintC3'
                    obj.computeMaxComplianceCthree();
                    obj.saveLoadCaseCthree();
                case 'complianceConstraintC4'
                    obj.computeMaxComlianceCfour();
                    obj.saveLoadCasefour();
                otherwise
                    
            end
        end

        function computeMaxComplianceCone(obj)
            run("bridgeLoadCasesInformation.m");
            obj.physicalProblem.boundaryConditions.changeBoundaryConditions(g.four',g.fourv');
            obj.physicalProblem.setC(obj.homogenizedVariablesComputer.C);
            obj.physicalProblem.computeStiffnessMatrix();
            obj.physicalProblem.solve();
            obj.targetConstraint = obj.computeInitialCompliance();
        end

        function computeMaxComplianceCtwo(obj)
            run("bridgeLoadCasesInformation.m");
            c1 = obj.computeCompliance(g.zero,g.zerov);
            c2 = obj.computeCompliance(g.four,g.fourv);
            c3 = obj.computeCompliance(g.eight,g.eightv);
            obj.targetConstraint = max([c1,c2,c3]);
        end

        function computeMaxComplianceCthree(obj)
            run("bridgeLoadCasesInformation.m");
            c1 = obj.computeCompliance(g.zero,g.zerov);
            c2 = obj.computeCompliance(g.one,g.onev);
            c3 = obj.computeCompliance(g.two,g.twov);
            c4 = obj.computeCompliance(g.three,g.threev);
            c5 = obj.computeCompliance(g.four,g.fourv);
            c6 = obj.computeCompliance(g.five,g.fivev);
            c7 = obj.computeCompliance(g.six,g.sixv);
            c8 = obj.computeCompliance(g.seven,g.sevenv);
            c9 = obj.computeCompliance(g.eight,g.eightv);
            obj.targetConstraint = max([c1,c2,c3,c4,c5,c6,c7,c8,c9]);
        end

        function computeMaxComplianceCfour(obj)
            run("bridgeLoadCasesInformation.m");
            c1 = obj.computeCompliance(obj,g.zero,g.zerov);
            c2 = obj.computeCompliance(obj,g.one,g.onev);
            c3 = obj.computeCompliance(obj,g.two,g.twov);
            c4 = obj.computeCompliance(obj,g.three,g.threev);
            obj.targetConstraint = max(c1,c2,c3,c4);
        end

        function saveLoadCaseCone(obj)
            run("bridgeLoadCasesInformation.m");
            obj.loadCase.neumann       = g.four;
            obj.loadCase.neumann_values = g.fourv;
        end

        function saveLoadCaseCtwo(obj)
            load = obj.loadCase.number;
            run("bridgeLoadCasesInformation.m");
            if load == 1
                obj.loadCase.neumann       = g.zero;
                obj.loadCase.neumann_values = g.zerov;
            elseif load == 2
                obj.loadCase.neumann = g.four;
                obj.loadCase.neumann_values = g.fourv;
            elseif load == 3
                obj.loadCase.neumann = g.eight;
                obj.loadCase.neumann_values = g.eightv;
            end
        end

        function saveLoadCaseCthree(obj)
            l = obj.loadCase.number;
            run("bridgeLoadCasesInformation.m");
            if l == 1
                obj.loadCase.neumann = g.zero;
                obj.loadCase.neumann_values = g.zerov;
            elseif l == 2
                obj.loadCase.neumann = g.one;
                obj.loadCase.neumann_values = g.onev;
            elseif l == 3
                obj.loadCase.neumann = g.two;
                obj.loadCase.neumann_values = g.twov;
            elseif l == 4
                obj.loadCase.neumann = g.three;
                obj.loadCase.neumann_values = g.threev;
            elseif l == 5
                obj.loadCase.neumann = g.four;
                obj.loadCase.neumann_values = g.fourv;
            elseif l == 6
                obj.loadCase.neumann = g.five;
                obj.loadCase.neumann_values = g.fivev;
            elseif l == 7
                obj.loadCase.neumann = g.six;
                obj.loadCase.neumann_values = g.sixv;
            elseif l == 8
                obj.loadCase.neumann = g.seven;
                obj.loadCase.neumann_values = g.sevenv;
            elseif l == 9
                obj.loadCase.neumann = g.eight;
                obj.loadCase.neumann_values = g.eightv;
            end
        end

        function saveLoadCaseCfour(obj)
            l = obj.loadCase.number;
            run("bridgeLoadCasesInformation.m");
            if l == 1
                obj.loadCase.neumann = g.zero;
                obj.loadCase.neumann_values = g.zerov;
            elseif l == 2
                obj.loadCase.neumann = g.one;
                obj.loadCase.neumann_values = g.onev;
            elseif l == 3
                obj.loadCase.neumann = g.two;
                obj.loadCase.neumann_values = g.twov;
            elseif l == 4
                obj.loadCase.neumann = g.three;
                obj.loadCase.neumann_values = g.threev;
            end
        end

        function c = computeInitialCompliance(obj)
            phy    = obj.physicalProblem;
            dvolum = phy.getDvolume()';
            stress = phy.variables.stress;
            strain = phy.variables.strain;
            ngaus  = size(strain,1);
            nelem  = size(strain,3);
            c      = zeros(nelem,ngaus);
            for igaus = 1:ngaus
                stressG    = squeeze(stress(igaus,:,:));
                strainG    = squeeze(strain(igaus,:,:));
                e          = stressG.*strainG;
                c(:,igaus) = c(:,igaus) + sum(e)';
            end
            obj.compliance = c;
            int            = c.*dvolum;
            c              = sum(int(:));
        end
        
    end

    methods (Access = protected)
        
        function solveState(obj)
            if obj.isFirstIter
                obj.computeTargetConstraint();
                obj.isFirstIter = false;
            end
            n  = obj.loadCase.neumann;
            nv = obj.loadCase.neumann_values;
            obj.physicalProblem.boundaryConditions.changeBoundaryConditions(n,nv)
            obj.physicalProblem.setC(obj.homogenizedVariablesComputer.C);
            obj.physicalProblem.computeStiffnessMatrix();
            obj.physicalProblem.solve();
        end
        
        function solveAdjoint(obj)
            obj.adjointProblem = obj.physicalProblem;
        end
        
        function computeFunctionValue(obj)
            phy    = obj.physicalProblem;
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
            if isempty(obj.value)
                obj.value = sum(int(:));
            else
                obj.value = sum(int(:));
                obj.value = obj.value - 0.7*obj.targetConstraint;
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