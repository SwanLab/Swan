classdef SettingsShapeFunctionalFactory < handle
    
    methods (Access = public)
        
        function s = create(obj,cParams)
            switch cParams.type
                case {'compliance','perimeter','volume','volumeConstraint',...
                        'chomog_CC','enforceCh_CCstar_L2','nonadjoint_compliance'}
                    s = SettingsShapeFunctional();
                case 'perimeterConstraint'
                    s = SettingsShFunc_PerimeterConstraint();
                    s.Perimeter_target = cParams.PerimeterTarget;
                case {'chomog_alphabeta','chomog_fraction'}
                    s = SettingsShFunc_Chomog();
                    s.alpha = cParams.alpha;
                    s.beta = cParams.beta;
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

