classdef Cost < handle
    properties
        ShapeFuncs
        nSF
        weights
        value
        gradient
        target_parameters
    end
    methods
        function obj = Cost(settings,weights)
            if isempty(weights)
                weights = ones(1,length(settings.cost));
            end
            obj.weights = weights;
            obj.nSF = length(settings.cost);
            for iSF = 1:obj.nSF
                settings_this = settings;
                %% !! PENDING TO BE IMPLEMENTED A FILTER PER EACH COST
                %                 settings_this.filter = settings.filter{iSF};
                switch settings.cost{iSF}
                    case 'compliance'
                        obj.ShapeFuncs{iSF} = ShFunc_Compliance(settings_this);
                    case 'perimeter'
                        obj.ShapeFuncs{iSF} = ShFunc_Perimeter(settings_this);
                    case 'chomog_alphabeta'
                        obj.ShapeFuncs{iSF} = ShFunc_Chomog_alphabeta(settings_this);
                    case 'chomog_fraction'
                        obj.ShapeFuncs{iSF} = ShFunc_Chomog_fraction(settings_this);
                    case 'nonadjoint_compliance'
                        obj.ShapeFuncs{iSF} = ShFunc_NonSelfAdjoint_Compliance(settings_this);
                    otherwise
                        error('Wrong cost name or not added to Cost Object')
                end
            end
        end
        
        function preProcess(obj,params)
            for iSF = 1:obj.nSF
                obj.ShapeFuncs{iSF}.filter.preProcess(params);
            end
        end
        
        function computef(obj, x, physicalProblem, interpolation)
            obj.value = 0;
            obj.gradient = zeros(length(x),1);
            for iSF = 1:length(obj.ShapeFuncs)
                obj.ShapeFuncs{iSF}.target_parameters = obj.target_parameters;
                obj.ShapeFuncs{iSF}.computef(x, physicalProblem, interpolation);
                obj.value = obj.value+obj.weights(iSF)*obj.ShapeFuncs{iSF}.value;
                obj.gradient = obj.gradient+obj.weights(iSF)*obj.ShapeFuncs{iSF}.gradient;
            end
        end
    end
end

