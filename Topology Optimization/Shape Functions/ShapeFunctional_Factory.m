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
            
            cParams.mesh = cParams.designVariable.mesh.innerMeshOLD;     
            cParams.filterParams.mesh = cParams.designVariable.mesh.innerMeshOLD;            
            cParams.filterParams.designVarType = cParams.designVariable.type;
            
            switch cParams.type
                case 'compliance'
                    sF = ShFunc_Compliance(cParams);
                case 'stressNorm'
                    sF = ShFunc_StressNorm(cParams);
                    %sF = ShFunc_StressNorm2(cParams);
                    %sF = ShFunc_StressNorm3(cParams);                    
                case 'perimeter'
                    cParams.filterParams.femSettings.isRobinTermAdded = true;
                    %cParams.designVariable = cParams.designVariable.value;
                    sF = ShFunc_Perimeter(cParams);
                case 'perimeterInterior'
                    cParams.filterParams.femSettings.isRobinTermAdded = false;
                    %cParams.designVariable = cParams.designVariable.value;                    
                    sF = ShFunc_Perimeter(cParams);
                case 'perimeterConstraint'
                    cParams.filterParams.femSettings.isRobinTermAdded = true;   
                    %cParams.designVariable = cParams.designVariable.value;                    
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
    
end

