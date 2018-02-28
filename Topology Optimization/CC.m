classdef CC < handle
    %CC Common core of Cost and Contraints subclasses
    
    properties
        ShapeFuncs
        nSF
        value
        gradient
        target_parameters
    end
    
    methods
        function obj = CC(settings_this,SF_list)
            obj.target_parameters = settings_this.target_parameters;
            obj.nSF = length(SF_list);
            for iSF = 1:obj.nSF
                
                %% !! PENDING TO BE DEFINED A FILTER FOR EACH SF
                % settings_this.filter = settings.filter{iSF};
                
                switch SF_list{iSF}
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
                    case 'volume'
                        obj.ShapeFuncs{iSF} = ShFunc_Volume(settings_this);
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
        
        function computef(obj, x, interpolation)
            obj.value = 0;
            obj.gradient = zeros(length(x),1);
            for iSF = 1:length(obj.ShapeFuncs)
                obj.ShapeFuncs{iSF}.target_parameters = obj.target_parameters;
                obj.ShapeFuncs{iSF}.computef(x,interpolation);
                obj.updateFields(iSF);
            end
        end
    end
end

