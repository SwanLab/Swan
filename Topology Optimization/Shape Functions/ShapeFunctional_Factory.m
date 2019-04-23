classdef ShapeFunctional_Factory
    
    properties (Access = private)
        settings
        cParams
        designVar
    end
    
    methods (Access = public)
        
        function shapeFunction = create(obj,type,settings,designVar)
            obj.settings = settings;
            obj.designVar = designVar;
            switch type
                case 'compliance'
                    obj.cParams = SettingsShapeFunctional();
                    obj.addParamsFromSettings();
                    shapeFunction = ShFunc_Compliance(obj.cParams);
                case 'perimeter'
                    obj.cParams = SettingsShapeFunctional();
                    obj.addParamsFromSettings();
                    shapeFunction = ShFunc_Perimeter(obj.cParams);
                case 'perimeterConstraint'
                    obj.cParams = SettingsShFunc_PerimeterConstraint();
                    obj.addParamsFromSettings();
                    obj.cParams.Perimeter_target = settings.Perimeter_target;
                    shapeFunction = Perimeter_constraint(obj.cParams);
                case 'chomog_alphabeta'
                    obj.cParams = SettingsShFunc_Chomog();
                    obj.addParamsFromSettings();
                    obj.cParams.alpha=settings.micro.alpha;
                    obj.cParams.beta=settings.micro.beta;
                    shapeFunction = ShFunc_Chomog_alphabeta(obj.cParams);
                case 'chomog_fraction'
                    obj.cParams = SettingsShFunc_Chomog();
                    obj.addParamsFromSettings();
                    obj.cParams.alpha = settings.micro.alpha;
                    obj.cParams.beta  = settings.micro.beta;
                    shapeFunction = ShFunc_Chomog_fraction(obj.cParams);
                case 'chomog_CC'
                    obj.cParams = SettingsShapeFunctional();
                    obj.addNamePtype()
                    shapeFunction = ShFunc_Chomog_CC();
                case 'enforceCh_CCstar_inf'
                    obj.cParams = SettingsShapeFunctional();
                    obj.addNamePtype()
                    for i=1:6
                        EnforceCh=ShFunc_Chomog_EnforceCh_CCstar_inf(settings,i);
                        if isequal(i,5) || isequal(i,4)
                            EnforceCh.setEpsilon(0);
                        end
                        shapeFunction=EnforceCh;
                        iSF = iSF+1;
                    end
                case 'enforceCh_CCstar_eq'
                    obj.cParams = SettingsShapeFunctional();
                    obj.addNamePtype()
                    obj.createFilterParams();
                    for i=1:6
                        shapeFunction=ShFunc_Chomog_EnforceCh_CCstar_eq(settings,i);
                        iSF = iSF+1;
                    end
                case 'enforceCh_CCstar_L2'
                    obj.cParams = SettingsShapeFunctional();
                    obj.addParamsFromSettings();
                    shapeFunction = ShFunc_Chomog_EnforceCh_CCstar_L2(obj.cParams);
                case 'nonadjoint_compliance'
                    obj.cParams = SettingsShapeFunctional();
                    obj.addParamsFromSettings();
                    shapeFunction = ShFunc_NonSelfAdjoint_Compliance(obj.cParams);
                case 'volume'
                    obj.cParams = SettingsShapeFunctional();
                    obj.addParamsFromSettings();
                    shapeFunction = ShFunc_Volume(obj.cParams);
                case 'volumeConstraint'
                    obj.cParams = SettingsShapeFunctional();
                    obj.addParamsFromSettings()
                    shapeFunction = Volume_constraint(obj.cParams);
                otherwise
                    error('Wrong cost name or not added to Cost Object')
            end
        end
        
    end
    
    methods (Access = private)
        
        function addParamsFromSettings(obj)
            obj.addNamePtype()
            obj.createMaterialInterpolationParams();
            obj.createFilterParams();
        end
        
        function addNamePtype(obj)
            obj.cParams.filename     = obj.settings.filename;
            obj.cParams.domainType   = obj.settings.ptype;
        end
        
        function createMaterialInterpolationParams(obj)
            %s = SettingsHomogenizedVarComputerFromInterpolation();
            s.type = obj.settings.homegenizedVariablesComputer;
            s.interpolation          = obj.settings.materialInterpolation;
            s.dim                    = obj.settings.pdim;
            s.typeOfMaterial         = obj.settings.material;
            s.constitutiveProperties = obj.settings.TOL;
            s.vademecumFileName      = obj.settings.vademecumFileName;
            obj.cParams.materialInterpolationParams = s;
        end
        
        function createFilterParams(obj)
            s = SettingsFilter();
            s.filterType = obj.settings.filter;
            s.designVar  = obj.designVar;
            obj.cParams.filterParams = s;
        end
        
    end
    
end

