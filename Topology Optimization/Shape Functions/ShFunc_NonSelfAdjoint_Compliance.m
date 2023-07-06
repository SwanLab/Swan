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

        function t = getTitlesToPlot(obj)
            t{1} = 'Compliance non scaled';
        end

        function fP = createPrintVariables(obj)
            fP{1}.type  = 'Elasticity';
            fP{2}.type  = 'Elasticity';
            fP{1}.name  = 'Primal';
            fP{2}.name  = 'Dual';
        end

        function v = getVariablesToPlot(obj)
            v{1} = obj.value*obj.value0;
        end

         function [fun, funNames] = getFunsToPlot(obj)
            mesh = obj.designVariable.mesh;
            phy = obj.physicalProblem;
            strain = phy.strainFun{1}; % !!!
            stress = phy.stressFun{1}; % !!!
            displ  = phy.uFun{1}; % !!!
            compl  = obj.value/obj.value0;

            quad = Quadrature.set(mesh.type);
            quad.computeQuadrature('CONSTANT');

            aa.mesh       = mesh;
            aa.quadrature = quad;
            aa.fValues    = permute(compl, [3 2 1]);
            complFun = FGaussDiscontinuousFunction(aa);

            fun      = {complFun, strain, stress, displ};
            funNames = {'NeumannDisplacementNRG', 'strain', 'stress', 'u'};
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
            obj.physicalProblem.solve();
        end

        function computeGradientValue(obj)
            eu    = obj.physicalProblem.variables.strain;
            ep    = obj.adjointProblem.variables.strain;
            nelem = size(eu,3);
            ngaus = size(eu,1);
            nstre = size(eu,2);
            g = zeros(nelem,ngaus,obj.nVariables);
            for igaus = 1:ngaus
                for istre = 1:nstre
                    for jstre = 1:nstre
                        eu_i = squeeze(eu(igaus,istre,:));
                        ep_j = squeeze(ep(igaus,jstre,:));
                        for ivar = 1:obj.nVariables
                            dCij = squeeze(obj.homogenizedVariablesComputer.dC(istre,jstre,ivar,:,igaus));
                            g(:,igaus,ivar) = g(:,igaus,ivar) - eu_i.*dCij.*ep_j;
                        end
                    end
                end
            end
            obj.gradient = g;
        end
        
        function solveAdjoint(obj)
            obj.adjointProblem.setC(obj.homogenizedVariablesComputer.C);
            obj.adjointProblem.solve();
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

    end

    methods (Access = private)

        function createAdjointProblem(obj,fileName)
            fAdj               = Preprocess.getBC_adjoint(fileName);
            a.fileName         = fileName;
            s                  = FemDataContainer(a);
            s.bc.pointload     = fAdj;
            obj.adjointProblem = FEM.create(s);
        end
        
        function f = computeDisplacementWeight(obj)
            nnode = obj.designVariable.mesh.nnodes;
            ndim  = obj.designVariable.mesh.ndim;
            ndof  = nnode*ndim;
            f = zeros(ndof,1);
            BC = obj.adjointProblem.boundaryConditions;
            if ~isempty(BC.neumann)
                f(obj.adjointProblem.boundaryConditions.neumann) = BC.neumann_values;
            end
        end
    end
end
