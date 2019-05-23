classdef SettingsShapeFunctionalFactory < handle
    
    methods (Access = public)
        
        function s = create(obj,cParams)
            switch cParams.type
                case {'compliance','perimeter','volume','volumeConstraint',...
                        'chomog_CC','enforceCh_CCstar_L2','nonadjoint_compliance','stressNorm'}
                    s = SettingsShapeFunctional();
                case 'perimeterConstraint'
                    s = SettingsShFunc_PerimeterConstraint();
                    s.Perimeter_target = cParams.PerimeterTarget;
                case {'chomog_alphabeta','chomog_fraction'}
                    s = SettingsShFunc_Chomog();
                    s.alpha = cParams.alpha;
                    s.beta = cParams.beta;
                case 'enforceCh_CCstar_inf'
                    error('Settings still not implemented');
                case 'enforceCh_CCstar_eq'
                    error('Settings still not implemented');
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

