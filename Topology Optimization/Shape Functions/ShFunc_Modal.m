classdef ShFunc_Modal < ShFunWithElasticPdes

    properties (Access = private)
        compliance
        fieldToPrint
        adjointProblem

    end

    properties (Access = public)
        mesh
        Ep1
        ePrint
        v
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
            obj.physicalProblem.solve(x); % x
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
            ep    = phy.variables.strain;
            derM = phy.derMass;
            obj.mesh = phy.mesh;
            ngaus  = size(ep,1);
            nstre  = size(ep,2);
            nelem  = size(ep,3);
            gK = zeros(nelem,ngaus,obj.nVariables);
            for igaus = 1:ngaus
                for istre = 1:nstre
                    for jstre = 1:nstre
                        eu_i = squeeze(ep(igaus,istre,:));
                        ep_j = squeeze(ep(igaus,jstre,:));
                        for ivar = 1:obj.nVariables
                            dCij = squeeze(obj.homogenizedVariablesComputer.dC(istre,jstre,ivar,:,igaus));
                            gK(:,igaus,ivar) = gK(:,igaus,ivar) + (-eu_i.*dCij.*ep_j);
                        end
                    end
                end
            end
            eigValues =  phy.variables.eigenValues;
            val = eigValues(1);
            obj.v = val*derM;
            % x = obj.regDesignVariable;
            g = gK - val*1;%derM; % x{1}
            obj.gradient = g;
%             ep1 = zeros(3,size(ep,3));
% 
%             for i=1 : obj.mesh.nelem
%             ep1(1,i) = ep(1,1,i);
%             ep1(2,i) = ep(1,2,i);
%             ep1(3,i) = ep(1,3,i);
%             end
%             obj.Ep1 = ep1; 
        end

        function f = getPdesVariablesToPrint(obj)
            f{1} = obj.getPdeVariableToPrint(obj.physicalProblem);
        end

    end

end

