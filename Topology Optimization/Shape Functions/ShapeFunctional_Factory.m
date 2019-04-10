classdef ShapeFunctional_Factory
    
    methods (Static, Access = public)
        
        function shapeFunction = create(type,settings,designVar)
            
            %             if ~isempty(settings.shFuncParamsName)
            %                 new_settings=SettingsShapeFunctional(settings.shFuncParamsName);
            %             else
            filterParams = SettingsFilter();
            filterParams.filterType = settings.filter;
            filterParams.designVar = designVar;
            
            materialInterpolationParams = SettingsInterpolation();
            materialInterpolationParams.interpolation=settings.method;
            materialInterpolationParams.dim=settings.pdim;
            materialInterpolationParams.typeOfMaterial=settings.material;
            materialInterpolationParams.constitutiveProperties=settings.TOL;
            
            d = SettingsShapeFunctional();
            
            d.filterParams = filterParams;
            d.materialInterpolationParams = materialInterpolationParams;
            d.filename = settings.filename;
            d.domainType = settings.ptype;
            
            new_settings=d;
            %             end
            
            
            switch type
                case 'compliance'
                    shapeFunction = ShFunc_Compliance(new_settings);
                case 'perimeter'
                    shapeFunction = ShFunc_Perimeter(new_settings);
                case 'perimeterConstraint'
                    d=SettingsShFunc_PerimeterConstraint();
                    d.filename=settings.filename;
                    d.domainType=settings.ptype;
                    d.filterParams.filterType=settings.filter;
                    d.filterParams.designVar = designVar;
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
                    d.filterParams.filterType=settings.filter;
                    d.filterParams.designVar = designVar;
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
                    d.filterParams.filterType=settings.filter;
                    d.filterParams.designVar = designVar;
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

