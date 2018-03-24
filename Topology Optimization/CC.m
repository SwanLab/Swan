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
                    case 'chomog_CC'
                        obj.ShapeFuncs{iSF} = ShFunc_Chomog_CC(settings_this);
                    case 'enforceCh_CCstar_inf'
                        for i=1:6
                            EnforceCh=ShFunc_Chomog_EnforceCh_CCstar_inf(settings_this,i);
                            if isequal(i,5) || isequal(i,4)
                                EnforceCh.setEpsilon(0);
                            end
                            obj.ShapeFuncs{iSF}=EnforceCh;
                            iSF = iSF+1;
                        end
                    case 'enforceCh_CCstar_eq'
                        for i=1:6
                            obj.ShapeFuncs{iSF}=ShFunc_Chomog_EnforceCh_CCstar_eq(settings_this,i);
                            iSF = iSF+1;
                        end
                    case 'enforceCh_CCstar_L2'
                        obj.ShapeFuncs{iSF}=ShFunc_Chomog_EnforceCh_CCstar_L2(settings_this);
                    case 'nonadjoint_compliance'
                        obj.ShapeFuncs{iSF} = ShFunc_NonSelfAdjoint_Compliance(settings_this);
                    case 'volume'
                        obj.ShapeFuncs{iSF} = ShFunc_Volume(settings_this);
                    otherwise
                        error('Wrong cost name or not added to Cost Object')
                end
            end
            obj.nSF = length(obj.ShapeFuncs);
        end
        
        function preProcess(obj,params)
            for iSF = 1:obj.nSF
                obj.ShapeFuncs{iSF}.filter.preProcess(params);
            end
        end
        
        function computef(obj, x)
            obj.value = 0;
            obj.gradient = zeros(length(x),1);
            for iSF = 1:length(obj.ShapeFuncs)
                obj.updateTargetParameters(iSF);
                obj.ShapeFuncs{iSF}.computef(x);
                obj.updateFields(iSF);
            end
        end
        
        function updateTargetParameters(obj,iSF)
            obj.ShapeFuncs{iSF}.target_parameters = obj.target_parameters;
            if isprop(obj.ShapeFuncs{iSF}.filter,'epsilon')
                obj.ShapeFuncs{iSF}.filter.epsilon=obj.target_parameters.epsilon;
            end
        end
    end
end

