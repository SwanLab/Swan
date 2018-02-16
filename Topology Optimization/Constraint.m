classdef Constraint < Shape_Functional
    properties 
        ShapeFuncs
        lambda
    end
    methods 
        function obj=Constraint(settings)
            for ifunc=1:length(settings.constraint)
                switch settings.constraint{ifunc}
                    case 'compliance'
                        obj.ShapeFuncs{ifunc}=ShFunc_Compliance(settings);
                    case 'perimeter'
                        obj.ShapeFuncs{ifunc}=ShFunc_Perimeter(settings);
                    case 'volume'
                        obj.ShapeFuncs{ifunc}=ShFunc_Volume(settings);
                    otherwise
                        error('Wrong constraint name or not added to Constraint Object')
                end
            end
        end
        function computef(obj, x, physicalProblem, interpolation,filter)
            for ifunc=1:length(obj.ShapeFuncs)
                obj.ShapeFuncs{ifunc}.target_parameters=obj.target_parameters;
                obj.ShapeFuncs{ifunc}.computef(x, physicalProblem, interpolation,filter);
                obj.value(ifunc,1)=obj.ShapeFuncs{ifunc}.value;
                obj.gradient(:,ifunc)=obj.ShapeFuncs{ifunc}.gradient;
            end
        end
    end
end