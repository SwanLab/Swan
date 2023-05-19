classdef ShFun_ComplianceTarget < TargetComputer
    methods (Access = public)
        function createShFunc(obj,cParams)
            obj.target = 2;
            cParams.type  = 'compliance';
            f = ShapeFunctional_Factory();
            obj.ShapeFunction = f.create(cParams);
        end
    end
end