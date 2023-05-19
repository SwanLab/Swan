classdef ShFun_ComplianceTarget < TargetComputer
    methods (Access = public)
        function createShFunc(obj,cParams)
            cParams.filterParams.femSettings.eta  = 0.75;
            cParams.filterParams.femSettings.beta = 1;
            obj.target = 2;
            cParams.type  = 'compliance';
            f = ShapeFunctional_Factory();
            obj.ShapeFunction = f.create(cParams);
        end
    end
end