classdef SettingsShapeFunctionalFactory < handle

    methods (Access = public)

        function s = create(obj,cParams)
            switch cParams.type
                case {'compliance','complianceConstraint','perimeter','perimeterInterior','volume','volumeConstraint',...
                        'chomog_CC','nonadjoint_compliance','stressNorm','anisotropicPerimeter2D','anisotropicPerimeterInterior2D'}
                    s = SettingsShapeFunctional(cParams);
                case 'perimeterConstraint'
                    s = SettingsShFunc_PerimeterConstraint(cParams);
                case {'chomog_alphabeta','chomog_fraction','enforceCh_CCstar_L2'}
                    s = SettingsShFunc_Chomog(cParams);
                case 'enforceCh_CCstar_inf'
                    error('Settings still not implemented');
                case 'enforceCh_CCstar_eq'
                    error('Settings still not implemented');
                otherwise
                    error('Wrong cost name or not added to Cost Object')
            end
        end

    end

end
