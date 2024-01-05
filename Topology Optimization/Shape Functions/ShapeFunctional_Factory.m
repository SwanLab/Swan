classdef ShapeFunctional_Factory < handle
    
    properties (Access = private)
        cParams
        designVar
        homogVarComputer
        targetParameters
    end
    
    methods (Access = public)
        
        function sF = create(obj,cParams)
            obj.designVar        = cParams.designVariable;

            if isprop(cParams,'homogVarComputer')
                obj.homogVarComputer = cParams.homogVarComputer;
            end
            if isprop(cParams,'targetParameters')
                obj.targetParameters = cParams.targetParameters;
            end
            
            if ~isempty(cParams.designVariable.mesh)
%                 if isprop(cParams.designVariable.mesh,'innerMeshOLD')
                mOld = cParams.designVariable.mesh;
                cParams.mesh = mOld;
                cParams.filterParams.mesh          = mOld;
                cParams.filterParams.designVarType = cParams.designVariable.type;
%                 end
            end
            
            switch cParams.type
                case 'compliance'
                    sF = ShFunc_Compliance(cParams);
                case {'complianceConstraintC1','complianceConstraintC2','complianceConstraintC3',...
                        'complianceConstraintC4'}
                    sF = ShFunc_Compliance_constraint(cParams);
                case 'complianceConstraint'
                    sF = ShFunc_ComplianceComparison_constraint(cParams);
                case 'stressNorm'
                    sF = ShFunc_StressNorm(cParams);
                case {'perimeter','perimeterInterior','anisotropicPerimeter2D','anisotropicPerimeterInterior2D'}
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
                case 'firstEignValue_functional'
                    sF = ShFunc_FirstEigenValue(cParams);
                case 'doubleEig'
                    sF = Sh_doubleEig(cParams);
                case 'volumeColumn'
                    sF = Sh_volumeColumn(cParams);
                case 'firstEigTop'
                    sF = Sh_firstEigTop(cParams);
                case 'ComplianceConstraintThreeFieldRhoE'
                    cParams.filterParams.femSettings.eta  = 0.5+0.25;
                    cParams.filterParams.femSettings.beta = 1;
                    cParams.filterParams.femSettings.shFunType = 'compliance';
                    sF = ShFunc_ComplianceWithBoundTarget(cParams);
                case 'ComplianceConstraintThreeFieldRhoI'
                    cParams.filterParams.femSettings.eta  = 0.5;
                    cParams.filterParams.femSettings.beta = 1;
                    cParams.filterParams.femSettings.shFunType = 'compliance';
                    sF = ShFunc_ComplianceWithBoundTarget(cParams);
                case 'ComplianceConstraintThreeFieldRhoD'
                    cParams.filterParams.femSettings.eta  = 0.5-0.25;
                    cParams.filterParams.femSettings.beta = 1;
                    cParams.filterParams.femSettings.shFunType = 'compliance';
                    sF = ShFunc_ComplianceWithBoundTarget(cParams);
                case 'VolumeConstraintRhoD'
                    cParams.filterParams.femSettings.eta  = 0.5;
                    cParams.filterParams.femSettings.beta = 1;
                    sF = Volume_constraintWithBound(cParams);
                case 'LinearBoundFunction'
                    sF = LinearBoundFunction(cParams);
                otherwise
                    error('Wrong cost name or not added to Cost Object')
            end
        end

    end

end

