classdef ShapeFunctional_Factory
    methods (Static, Access = public)
        function shapeFunction = create(type,settings)
            filterParams=SettingsFilter();
            filterParams.optimizer=settings.optimizer;
            filterParams.filter=settings.filter;
            materialInterpolationParams = SettingsInterpolation();
            materialInterpolationParams.interpolation=settings.method;
            materialInterpolationParams.dim=settings.pdim;
            materialInterpolationParams.typeOfMaterial=settings.material;
            materialInterpolationParams.constitutiveProperties=settings.TOL;
            d=SettingsShapeFunctional(filterParams,materialInterpolationParams);
            d.filename=settings.filename;
            d.domainType=settings.ptype;

            new_settings=d;
            switch type
                case 'compliance'
                    shapeFunction = ShFunc_Compliance(new_settings);
                case 'perimeter'
                    shapeFunction = ShFunc_Perimeter(new_settings);
                case 'perimeterConstraint'
                    d=SettingsShFunc_PerimeterConstraint();
                    d.filename=settings.filename;
                    d.domainType=settings.ptype;
                    d.filterParams.optimizer=settings.optimizer;
                    d.filterParams.filter=settings.filter;
                    d.materialInterpolationParams.interpolation=settings.method;
                    d.materialInterpolationParams.dim=settings.pdim;
                    d.materialInterpolationParams.typeOfMaterial=settings.material;
                    d.materialInterpolationParams.constitutiveProperties=settings.TOL;
                    d.Perimeter_target=settings.Perimeter_target;
                    shapeFunction = Perimeter_constraint(d);
                case 'chomog_alphabeta'
                    d=SettingsShFunc_Chomog();
                    d.filename=settings.filename;
                    d.domainType=settings.ptype;
                    d.filterParams.optimizer=settings.optimizer;
                    d.filterParams.filter=settings.filter;
                    d.materialInterpolationParams.interpolation=settings.method;
                    d.materialInterpolationParams.dim=settings.pdim;
                    d.materialInterpolationParams.typeOfMaterial=settings.material;
                    d.materialInterpolationParams.constitutiveProperties=settings.TOL;
                    d.alpha=settings.micro.alpha;
                    d.beta=settings.micro.beta;
                    settings=d;
                    shapeFunction = ShFunc_Chomog_alphabeta(settings);
                case 'chomog_fraction'
                    d=SettingsShFunc_Chomog();
                    d.filename=settings.filename;
                    d.domainType=settings.ptype;
                    d.filterParams.optimizer=settings.optimizer;
                    d.filterParams.filter=settings.filter;
                    d.materialInterpolationParams.interpolation=settings.method;
                    d.materialInterpolationParams.dim=settings.pdim;
                    d.materialInterpolationParams.typeOfMaterial=settings.material;
                    d.materialInterpolationParams.constitutiveProperties=settings.TOL;
                    d.alpha=settings.micro.alpha;
                    d.beta=settings.micro.beta;
                    settings=d;
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
                    shapeFunction = ShFunc_NonSelfAdjoint_Compliance(new_settings);
                case 'volume'
                    shapeFunction = ShFunc_Volume(settings);
                case 'volumeConstraint'
                    shapeFunction = Volume_constraint(new_settings);
                otherwise
                    error('Wrong cost name or not added to Cost Object')
            end
        end
    end
end

