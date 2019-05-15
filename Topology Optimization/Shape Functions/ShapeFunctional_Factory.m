classdef ShapeFunctional_Factory < handle
    
    properties (Access = private)
        cParams
        designVar
        homogVarComputer
        targetParameters
    end
    
    methods (Access = public)
        
        function sF = create(obj,cParams)
            obj.designVar = cParams.designVariable;
            obj.homogVarComputer = cParams.homogVarComputer;
            obj.targetParameters = cParams.targetParameters;
            
            cParams.filterParams.designVar = cParams.designVariable;
            switch cParams.type
                case 'compliance'
                    sF = ShFunc_Compliance(cParams);
                case 'perimeter'
                    sF = ShFunc_Perimeter(cParams);
                case 'perimeterConstraint'
                    sF = Perimeter_constraint(cParams);
                case 'chomog_alphabeta'
                    sF = ShFunc_Chomog_alphabeta(cParams);
                case 'chomog_fraction'
                    sF = ShFunc_Chomog_fraction(cParams);
                case 'chomog_CC'
                    sF = ShFunc_Chomog_CC();
                case 'enforceCh_CCstar_inf'
                    for i = 1:6
                        EnforceCh = ShFunc_Chomog_EnforceCh_CCstar_inf(settings,i);
                        if isequal(i,5) || isequal(i,4)
                            EnforceCh.setEpsilon(0);
                        end
                        sF  = EnforceCh;
                        iSF = iSF+1;
                    end
                case 'enforceCh_CCstar_eq'
                    for i = 1:6
                        sF  = ShFunc_Chomog_EnforceCh_CCstar_eq(settings,i);
                        iSF = iSF+1;
                    end
                case 'enforceCh_CCstar_L2'
                    sF = ShFunc_Chomog_EnforceCh_CCstar_L2(cParams);
                case 'nonadjoint_compliance'
                    sF = ShFunc_NonSelfAdjoint_Compliance(cParams);
                case 'volume'
                    sF = ShFunc_Volume(cParams);
                case 'volumeConstraint'
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

