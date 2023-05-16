classdef ShFunc_MorphologyBasedCompliance < ShapeFunctional
    properties (Access = private)
        compplianceE
        compplianceI
        compplianceD
        factory
    end
    methods (Access = public)
        function obj = ShFunc_MorphologyBasedCompliance(cParams)
            obj.createCompliances(cParams);
        end
        function updateTargetParameters(obj)
             obj.compplianceE.updateTargetParameters();
             obj.compplianceD.updateTargetParameters();
             obj.compplianceI.updateTargetParameters();
        end 
        function computeFunctionAndGradient(obj)
             obj.compplianceE.computeFunctionAndGradient();
             obj.compplianceD.computeFunctionAndGradient();
             obj.compplianceI.computeFunctionAndGradient();
             obj.order();
        end 
        function computeFunction(obj)
            obj.compplianceE.computeFunction();
            obj.compplianceI.computeFunction();
            obj.compplianceD.computeFunction();
            obj.order();
        end
        function computeGradient(obj)
            obj.compplianceE.computeFunction();
            obj.compplianceI.computeFunction();
            obj.compplianceD.computeFunction();
            obj.order();
        end
        function f = getTitlesToPlot(obj)
            f = {'MorphologyBasedCompliance'};
        end
    end 
    methods (Access = protected)
%         function init(obj,cParams)
%             obj.createCompliances(cParams);
%         end
        function cParams = createCompliances(obj,cParams)
            cParams.filterParams.femSettings.beta = 1;
            obj.createComplianceE(cParams);
            obj.createComplianceI(cParams);
            obj.createComplianceD(cParams);
        end
        function createComplianceE(obj,cParams)
            cParams.type = 'ComplianceConstraintThreeFieldRhoE';
            cParams.filterParams.femSettings.eta  = 0.25;
            obj.compplianceE = obj.create(cParams);
        end
        function createComplianceI(obj,cParams)
            cParams.type = 'ComplianceConstraintThreeFieldRhoI';
            cParams.filterParams.femSettings.eta  = 0.5;
            obj.compplianceI = obj.create(cParams);
        end
        function createComplianceD(obj,cParams)
            cParams.type = 'ComplianceConstraintThreeFieldRhoD';
            cParams.filterParams.femSettings.eta  = 0.75;
            obj.compplianceD = obj.create(cParams);
        end
        function order(obj)
            obj.value = [obj.compplianceE.value;obj.compplianceI.value;obj.compplianceD.value];
            obj.gradient = [obj.compplianceE.gradient;obj.compplianceI.gradient;obj.compplianceD.gradient];
        end 
    end

end

