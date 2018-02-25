
classdef Constraint < handle
    properties
        ShapeFuncs
        weights
        value
        gradient
        target_parameters
        nSF
        lambda
    end
    methods
        function obj=Constraint(settings)
            obj.nSF = length(settings.constraint);
            for iSF=1:obj.nSF
                switch settings.constraint{iSF}
                    case 'compliance'
                        obj.ShapeFuncs{iSF}=ShFunc_Compliance(settings);
                    case 'perimeter'
                        obj.ShapeFuncs{iSF}=ShFunc_Perimeter(settings);
                    case 'volume'
                        obj.ShapeFuncs{iSF}=ShFunc_Volume(settings);
                    otherwise
                        error('Wrong constraint name or not added to Constraint Object')
                end
            end
        end
        function preProcess(obj,params)
            for iSF = 1:obj.nSF
                obj.ShapeFuncs{iSF}.filter.preProcess(params);
            end
        end
        function computef(obj, x, physicalProblem, interpolation)
            for iSF=1:length(obj.ShapeFuncs)
                obj.ShapeFuncs{iSF}.target_parameters=obj.target_parameters;
                obj.ShapeFuncs{iSF}.computef(x, physicalProblem, interpolation);
                obj.value(iSF,1)=obj.ShapeFuncs{iSF}.value;
                obj.gradient(:,iSF)=obj.ShapeFuncs{iSF}.gradient;
            end
        end
    end
end