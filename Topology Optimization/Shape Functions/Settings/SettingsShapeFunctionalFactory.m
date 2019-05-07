classdef SettingsShapeFunctionalFactory < handle
    
    methods (Access = public)
        
        function s = create(obj,cParams,settings)
            switch cParams.type
                case 'compliance'
                    s = SettingsShapeFunctional();
                case 'perimeter'
                    s = SettingsShapeFunctional();
                case 'perimeterConstraint'
                    s = SettingsShFunc_PerimeterConstraint();
                    s.Perimeter_target = settings.Perimeter_target;
                case 'chomog_alphabeta'
                    s = SettingsShFunc_Chomog();
                    s.alpha = settings.micro.alpha;
                    s.beta = settings.micro.beta;
                case 'chomog_fraction'
                    s = SettingsShFunc_Chomog();
                    s.alpha = settings.micro.alpha;
                    s.beta = settings.micro.beta;
                case 'chomog_CC'
                    s = SettingsShapeFunctional();
                case 'enforceCh_CCstar_inf'
                    % !! INCOMPLET !!
                    obj.cParams = SettingsShapeFunctional();
                    obj.addNamePtype()
                    for i = 1:6
                        EnforceCh = ShFunc_Chomog_EnforceCh_CCstar_inf(obj.settings,i);
                        if isequal(i,5) || isequal(i,4)
                            EnforceCh.setEpsilon(0);
                        end
                        s = EnforceCh;
                        iSF = iSF+1;
                    end
                case 'enforceCh_CCstar_eq'
                    % !! INCOMPLET !!
                    obj.cParams = SettingsShapeFunctional();
                    obj.addNamePtype()
                    obj.createFilterParams();
                    for i = 1:6
                        s = ShFunc_Chomog_EnforceCh_CCstar_eq(obj.settings,i);
                        iSF = iSF+1;
                    end
                case 'enforceCh_CCstar_L2'
                    s = SettingsShapeFunctional();
                case 'nonadjoint_compliance'
                    s = SettingsShapeFunctional();
                case 'volume'
                   s = SettingsShapeFunctional();
                case 'volumeConstraint'
                    s = SettingsShapeFunctional();
                otherwise
                    error('Wrong cost name or not added to Cost Object')
            end
            
            s = obj.setCommonParams(s,cParams);
            
        end
        
    end
    
    methods (Access = private)
        
        function s = setCommonParams(obj,s,cParams)
            s.filename = cParams.filename;
            s.scale = cParams.scale;
            s.filterParams = cParams.filterParams;
            s.type = cParams.type;
        end
        
    end
    
end

