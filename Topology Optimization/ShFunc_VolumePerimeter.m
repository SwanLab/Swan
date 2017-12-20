classdef ShFunc_VolumePerimeter< Shape_Functional
    properties 
        volume
        perimeter
    end
    methods 
        function obj=ShFunc_VolumePerimeter(settings)
        obj@Shape_Functional(settings);
        obj.volume=ShFunc_Volume(settings);
        obj.perimeter=ShFunc_Perimeter(settings);
        end
    end
    methods 
        function computef(obj, x, physicalProblem, interpolation,filter)
            obj.volume.target_parameters=obj.target_parameters;
            %obj.perimeter.target_parameters=obj.target_parameters;
            obj.volume.computef(x, physicalProblem, interpolation,filter);
            obj.perimeter.computef(x, physicalProblem, interpolation,filter);
            obj.value(1,1)=obj.volume.value;
            obj.value(2,1)=obj.perimeter.value;
            obj.gradient(:,1)=obj.volume.gradient;
            obj.gradient(:,2)=obj.perimeter.gradient;
        end
    end
end