classdef ShapeFunctional_Factory
    methods (Static, Access = public)
        function shapeFunction = create(type,settings)
            switch type
                case 'compliance'
                    shapeFunction = ShFunc_Compliance(settings);
                case 'perimeter'
                    shapeFunction = ShFunc_Perimeter(settings);
                case 'perimeterConstraint'
                    shapeFunction = Perimeter_constraint(settings);
                case 'chomog_alphabeta'
                    shapeFunction = ShFunc_Chomog_alphabeta(settings);
                case 'chomog_fraction'
                    shapeFunction = ShFunc_Chomog_fraction(settings);
                case 'chomog_CC'
                    shapeFunction = ShFunc_Chomog_CC(settings);
                case 'enforceCh_CCstar_inf'
                    for i=1:6
                        EnforceCh=ShFunc_Chomog_EnforceCh_CCstar_inf(settings,i);
                        if isequal(i,5) || isequal(i,4)
                            EnforceCh.setEpsilon(0);
                        end
                        shapeFunction=EnforceCh;
                        iSF = iSF+1;
                    end
                case 'enforceCh_CCstar_eq'
                    for i=1:6
                        shapeFunction=ShFunc_Chomog_EnforceCh_CCstar_eq(settings,i);
                        iSF = iSF+1;
                    end
                case 'enforceCh_CCstar_L2'
                    shapeFunction=ShFunc_Chomog_EnforceCh_CCstar_L2(settings);
                case 'nonadjoint_compliance'
                    shapeFunction = ShFunc_NonSelfAdjoint_Compliance(settings);
                case 'volume'
                    shapeFunction = ShFunc_Volume(settings);
                case 'volumeConstraint'
                    shapeFunction = Volume_constraint(settings);
                otherwise
                    error('Wrong cost name or not added to Cost Object')
            end
        end
    end
end

