classdef ShapeFunctional_Factory
    
    properties (Access = private)
        settings
        cParams
        designVar
        homogVarComputer
        targetParameters
    end
    
    methods (Access = public)
        
        function sF = create(obj,cParams,settings)
            obj.settings = settings;
            obj.designVar = cParams.designVariable;
            obj.homogVarComputer = cParams.homogVarComputer;
            obj.targetParameters = cParams.targetParameters;
            switch cParams.type
                case 'compliance'
                    obj.cParams = SettingsShapeFunctional();
                    obj.addParamsFromSettings();
                    sF = ShFunc_Compliance(obj.cParams);
                case 'perimeter'
                    obj.cParams = SettingsShapeFunctional();
                    obj.addParamsFromSettings();
                    sF = ShFunc_Perimeter(obj.cParams);
                case 'perimeterConstraint'
                    cParams.filterParams.designVar = cParams.designVariable;
                    sF = Perimeter_constraint(cParams);
                case 'chomog_alphabeta'
                    cParams.filterParams.designVar = cParams.designVariable;
                    sF = ShFunc_Chomog_alphabeta(cParams);
                case 'chomog_fraction'
                    obj.cParams = SettingsShFunc_Chomog();
                    obj.addParamsFromSettings();
                    obj.cParams.alpha = obj.settings.micro.alpha;
                    obj.cParams.beta  = obj.settings.micro.beta;
                    sF = ShFunc_Chomog_fraction(obj.cParams);
                case 'chomog_CC'
                    obj.cParams = SettingsShapeFunctional();
                    obj.addNamePtype()
                    sF = ShFunc_Chomog_CC();
                case 'enforceCh_CCstar_inf'
                    obj.cParams = SettingsShapeFunctional();
                    obj.addNamePtype()
                    for i=1:6
                        EnforceCh=ShFunc_Chomog_EnforceCh_CCstar_inf(obj.settings,i);
                        if isequal(i,5) || isequal(i,4)
                            EnforceCh.setEpsilon(0);
                        end
                        sF=EnforceCh;
                        iSF = iSF+1;
                    end
                case 'enforceCh_CCstar_eq'
                    obj.cParams = SettingsShapeFunctional();
                    obj.addNamePtype()
                    obj.createFilterParams();
                    for i=1:6
                        sF=ShFunc_Chomog_EnforceCh_CCstar_eq(obj.settings,i);
                        iSF = iSF+1;
                    end
                case 'enforceCh_CCstar_L2'
                    obj.cParams = SettingsShapeFunctional();
                    obj.addParamsFromSettings();
                    sF = ShFunc_Chomog_EnforceCh_CCstar_L2(obj.cParams);
                case 'nonadjoint_compliance'
                    obj.cParams = SettingsShapeFunctional();
                    obj.addParamsFromSettings();
                    sF = ShFunc_NonSelfAdjoint_Compliance(obj.cParams);
                case 'volume'
                    obj.cParams = SettingsShapeFunctional();
                    obj.addParamsFromSettings();
                    sF = ShFunc_Volume(obj.cParams);
                case 'volumeConstraint'
                    cParams.filterParams.designVar = cParams.designVariable;
                    sF = Volume_constraint(cParams);
                otherwise
                    error('Wrong cost name or not added to Cost Object')
            end
        end
        
    end
    
    methods (Access = private)
        
        function addParamsFromSettings(obj)
            obj.addNamePtype()
            obj.cParams.homogVarComputer = obj.homogVarComputer;
            obj.cParams.designVariable   = obj.designVar;
            obj.cParams.targetParameters = obj.targetParameters;
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

