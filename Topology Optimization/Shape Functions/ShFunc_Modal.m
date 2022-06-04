classdef ShFunc_Modal < ShFunWithElasticPdes
    
    properties (Access = private)
        compliance
        fieldToPrint
        adjointProblem
    end
    
    methods (Access = public)
        
        function obj = ShFunc_Modal(cParams)
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
            obj.physicalProblem.setC(obj.homogenizedVariablesComputer.C) % (:,:,7200,4); cmat
%             obj.physicalProblem.computeVariables();
            x = obj.regDesignVariable;
            obj.physicalProblem.solve(x);
        end
        
        function solveAdjoint(obj)
            obj.adjointProblem = obj.physicalProblem;
        end
        
        function computeFunctionValue(obj)
            phy = obj.physicalProblem;
            lambda = phy.variables.eigenValues;
            gamma = obj.designVariable.getFirstEigenMode();
            fx = gamma-lambda(1);
            obj.value = fx;
        end
        
        function computeGradientValue(obj)
             phy = obj.physicalProblem;
             mesh = phy.mesh;
             eigModes = phy.variables.eigenModes;
             eigValues =  phy.variables.eigenValues;
             mod = eigModes(:,1);
             val = eigValues(1);
             ngaus = 1; % pending to generalize
             nstre = 3; % pending to generalize
             nelem  = mesh.nelem;
             g_K = zeros(nelem,ngaus,obj.nVariables);
             for igaus = 1: ngaus
                 for istre = 1:nstre
                     for jstre = 1:nstre
                         for ivar = 1; obj.nVariables
                             dCij = squeeze(obj.homogenizedVariablesComputer.dC(istre,jstre,ivar,:,igaus));
                             g_K(:,igaus,ivar) = g_K(:,igaus,ivar) + (mod.*dCij.*mod);
                         end
                     end
                 end
                 g = (mod.*g_K.*mod - val*ones(nelem,1));
                obj.gradient = g;
             end


        end
        
        function f = getPdesVariablesToPrint(obj)
            f{1} = obj.getPdeVariableToPrint(obj.physicalProblem);
        end
        
    end
    
end

