classdef ShFunc_NonSelfAdjoint_Compliance < ShFunWithElasticPdes

    properties (Access = private)
        adjointProblem
    end

    methods (Access = public)

        function obj = ShFunc_NonSelfAdjoint_Compliance(cParams)
            cParams.filterParams.quadratureOrder = 'LINEAR';
            obj.init(cParams);
            fileName = cParams.femSettings.fileName;
            obj.createEquilibriumProblem(fileName);
            obj.createAdjointProblem(fileName);
            obj.createOrientationUpdater();                     
        end

    end

    methods (Access = protected)

        function computeFunctionValue(obj)
            u = obj.physicalProblem.variables.d_u;
            f = obj.computeDisplacementWeight();
            obj.value = f'*u;
        end
        
        function solveState(obj)
            obj.physicalProblem.setC(obj.homogenizedVariablesComputer.C);
            obj.physicalProblem.computeVariables();            
        end        

        function computeGradientValue(obj)
            eu    = obj.physicalProblem.variables.strain;
            ep    = obj.adjointProblem.variables.strain;
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
                            dCij = squeeze(obj.homogenizedVariablesComputer.dC(istre,jstre,ivar,:));
                            g(:,igaus,ivar) = eu_i.*dCij.*ep_j;
                        end
                    end
                end
            end
            obj.gradient = g;
        end
        
        function solveAdjoint(obj)
            obj.adjointProblem.setC(obj.homogenizedVariablesComputer.C);
            obj.adjointProblem.computeVariables();
        end
        
        function f = getPdesVariablesToPrint(obj)
            f{1} = obj.getPdeVariablesToPrint(obj.physicalProblem);
            f{2} = obj.getPdeVariablesToPrint(obj.adjointProblem);
        end        
        
        function fP = addPrintableVariables(obj)
            phy = obj.getPdesVariablesToPrint();
            fP{1}.value = phy{1};
            fP{2}.value = obj.compliance/obj.value0;
            fP{3}.value = obj.designVariable.alpha;
            fP{4}.value = abs(obj.designVariable.alpha);
            fP{5}.value = obj.getRegularizedDesignVariable();
            fP{6}.value = obj.homogenizedVariablesComputer.addPrintableVariables(obj.designVariable);
        end         
        
        function fP = createPrintVariables(obj)
            fP{1}.type  = 'Elasticity';
            fP{2}.type  = 'Elasticity';
            fP{1}.name  = 'Primal';            
            fP{2}.name  = 'Dual';            
        end        

    end

    methods (Access = private)

        function createAdjointProblem(obj,fileName)
            fAdj = Preprocess.getBC_adjoint(fileName);
            obj.adjointProblem = FEM.create(fileName);
            [dof,dofVal] = obj.adjointProblem.dof.get_dof_conditions(fAdj,obj.adjointProblem.dof.nunkn);
            obj.adjointProblem.dof.neumann = dof;
            obj.adjointProblem.dof.neumann_values = -dofVal;
        end
        
        function f = computeDisplacementWeight(obj)
            f = zeros(obj.adjointProblem.dof.ndof,1);
            if ~isempty(obj.adjointProblem.dof.neumann)
                f(obj.adjointProblem.dof.neumann) = obj.adjointProblem.dof.neumann_values;
            end
            f = -f;
        end
    end
end
