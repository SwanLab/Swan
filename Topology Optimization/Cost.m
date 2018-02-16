classdef Cost < Shape_Functional
    properties 
        ShapeFuncs
        multipliers
    end
    methods 
        function obj=Cost(settings,multipliers)
            if isempty(multipliers)
                multipliers=ones(1,length(settings.cost));
            end
            obj.multipliers=multipliers;
            for ifunc=1:length(settings.cost)
                switch settings.cost{ifunc}
                    case 'compliance'
                        obj.ShapeFuncs{ifunc}=ShFunc_Compliance(settings);
                    case 'perimeter'
                        obj.ShapeFuncs{ifunc}=ShFunc_Perimeter(settings);
                    case 'chomog_alphabeta'
                        obj.ShapeFuncs{ifunc}=ShFunc_Chomog_alphabeta(settings);
                    case 'chomog_fraction'
                        obj.ShapeFuncs{ifunc}=ShFunc_Chomog_fraction(settings);
                    otherwise
                        error('Wrong cost name or not added to Cost Object')
                end
            end
        end
        function computef(obj, x, physicalProblem, interpolation,filter)
            obj.value=0;
            obj.gradient=zeros(length(x),1);
            for ifunc=1:length(obj.ShapeFuncs)
                obj.ShapeFuncs{ifunc}.target_parameters=obj.target_parameters;
                obj.ShapeFuncs{ifunc}.computef(x, physicalProblem, interpolation,filter);
                obj.value=obj.value+obj.multipliers(ifunc)*obj.ShapeFuncs{ifunc}.value;
                obj.gradient=obj.gradient+obj.multipliers(ifunc)*obj.ShapeFuncs{ifunc}.gradient;
            end
        end
    end
end

    