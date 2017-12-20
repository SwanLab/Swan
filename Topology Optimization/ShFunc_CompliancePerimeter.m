classdef ShFunc_CompliancePerimeter< Shape_Functional
    properties 
        func
        compliance
        perimeter
        lambda
    end
    methods 
        function obj=ShFunc_CompliancePerimeter(settings)
        obj@Shape_Functional(settings);
        obj.compliance=ShFunc_Compliance(settings);
        obj.perimeter=ShFunc_Perimeter(settings);
        obj.lambda=settings.perimeter.lambda;
        end
    end
    methods 
        function computef(obj, x, physicalProblem, interpolation,filter)
            obj.compliance.computef(x, physicalProblem, interpolation,filter);
            obj.perimeter.computef(x, physicalProblem, interpolation,filter);
            obj.value=obj.compliance.value+obj.lambda*obj.perimeter.value;
            obj.gradient=obj.compliance.gradient+obj.lambda*obj.perimeter.gradient;
        end
    end
end