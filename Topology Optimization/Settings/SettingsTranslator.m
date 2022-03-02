classdef SettingsTranslator < handle

    properties (GetAccess = public, SetAccess = private)
        fileName
        filePath
    end

    properties (Access = private)
        oldSettings
        propList
        nProps
    end

    methods (Access = public)

        function translate(obj,oldSettings)
            obj.init(oldSettings);
            for i = 1:obj.nProps
                prop = obj.propList{i};
                value = oldSettings.(prop);
                if strcmp(prop,'filename')
                    s.problemData.femData.fileName = value;
                elseif strcmp(prop,'homegenizedVariablesComputer')
                    s.homogenizedVarComputerSettings.type = value;
                    s = obj.translateHomogenizedVariablesComputer(s);
                elseif strcmp(prop,'designVariable')
                    s.designVarSettings.type = value;
                elseif strcmp(prop,'initial_case')
                    s.designVarSettings.initialCase = value;
                elseif strcmp(prop,'m1')
                    if ~isempty(value)
                        s.designVarSettings.creatorSettings.m1 = value;
                    end
                elseif strcmp(prop,'m2')
                    if ~isempty(value)
                        s.designVarSettings.creatorSettings.m2 = value;
                    end
                elseif strcmp(prop,'alpha0')
                    if ~isempty(value)
                        s.designVarSettings.creatorSettings.alpha0 = value;
                    end
                elseif strcmp(prop,'rho0')
                    if ~isempty(value)
                        s.designVarSettings.creatorSettings.rho = value;
                    end
                elseif strcmp(prop,'levelSetDataBase')
                    s = obj.translateLevelSetCreator(s);
                elseif strcmp(prop,'nsteps')
                    s.incrementalSchemeSettings.nSteps = value;
                elseif strcmp(prop,'printIncrementalIter')
                    s.incrementalSchemeSettings.shallPrintIncremental = value;
                elseif strcmp(prop,'Vfrac_initial')
                    s.incrementalSchemeSettings.targetParamsSettings.VfracInitial = value;
                elseif strcmp(prop,'Vfrac_final')
                    s.incrementalSchemeSettings.targetParamsSettings.VfracFinal = value;
                elseif strcmp(prop,'optimality_initial')
                    s.incrementalSchemeSettings.targetParamsSettings.optimalityInitial = value;
                elseif strcmp(prop,'optimality_final')
                    s.incrementalSchemeSettings.targetParamsSettings.optimalityFinal = value;
                elseif strcmp(prop,'constr_initial')
                    s.incrementalSchemeSettings.targetParamsSettings.constrInitial = value;
                elseif strcmp(prop,'constr_final')
                    s.incrementalSchemeSettings.targetParamsSettings.constrFinal = value;
                elseif strcmp(prop,'stressNormExponent_initial')
                    s.incrementalSchemeSettings.targetParamsSettings.stressNormExponentInitial = value;
                elseif strcmp(prop,'stressNormExponent_final')
                    s.incrementalSchemeSettings.targetParamsSettings.stressNormExponentFinal = value;
                elseif strcmp(prop,'cost')
                    for k = 1:length(value)
                        s.costSettings.shapeFuncSettings{k}.type = value{k};
                        if strcmp(value{k},'chomog_alphabeta') || strcmp(value{k},'chomog_fraction')
                            if isprop(oldSettings,'micro')
                                s.costSettings.shapeFuncSettings{k}.alpha = oldSettings.micro.alpha;
                                s.costSettings.shapeFuncSettings{k}.beta = oldSettings.micro.beta;
                            end
                        elseif strcmp(value{k},'enforceCh_CCstar_L2')
                            s.costSettings.shapeFuncSettings{k}.ChTarget.type = oldSettings.selectiveC_Cstar;

                        elseif strcmp(value{k},'perimeterConstraint')
                            if isprop(oldSettings,'Perimeter_target')
                                s.costSettings.shapeFuncSettings{k}.perimeterTarget = oldSettings.Perimeter_target;
                            end
                        end
                    end
                elseif strcmp(prop,'weights')
                    s.costSettings.weights = value;
                elseif strcmp(prop,'filter')
                    s.costSettings.filterType = value;
                    s.constraintSettings.filterType = value;
                elseif strcmp(prop,'constraint')
                    for k = 1:length(value)
                        s.constraintSettings.shapeFuncSettings{k}.type = value{k};
                        if strcmp(value{k},'perimeterConstraint')
                            if isprop(oldSettings,'Perimeter_target')
                                s.constraintSettings.shapeFuncSettings{k}.perimeterTarget = oldSettings.Perimeter_target;
                            end
                        end
                    end
                elseif strcmp(prop,'constraintDomainNotOptimizable')
                    for k = 1:numel(value)
                        s.constraintSettings.shapeFuncSettings{k}.domainNotOptimizable = value{k};
                    end
                elseif strcmp(prop,'costDomainNotOptimizable')
                    for k = 1:numel(value)
                        s.costSettings.shapeFuncSettings{k}.domainNotOptimizable = value{k};
                    end
                elseif strcmp(prop,'optimizer')
                    s.optimizerSettings.type = value;
                elseif strcmp(prop,'constraint_case')
                    s.optimizerSettings.constraintCase = value;
                elseif strcmp(prop,'optimizerUnconstrained')
                    s.optimizerSettings.uncOptimizerSettings.type = value;
                elseif strcmp(prop,'e2')
                    s.optimizerSettings.uncOptimizerSettings.e2 = value;
                elseif strcmp(prop,'ub')
                    s.optimizerSettings.uncOptimizerSettings.ub = value;
                elseif strcmp(prop,'lb')
                    s.optimizerSettings.uncOptimizerSettings.lb = value;
                elseif strcmp(prop,'line_search_initiator')
                    s.optimizerSettings.uncOptimizerSettings.lineSearchSettings.lineSearchInitiatorSettings.type = value;
                elseif strcmp(prop,'HJiter0')
                    s.optimizerSettings.uncOptimizerSettings.lineSearchSettings.HJiter0 = value;
                elseif strcmp(prop,'incrementFactor')
                    s.optimizerSettings.uncOptimizerSettings.lineSearchSettings.lineSearchInitiatorSettings.incrementFactor = value;
                elseif strcmp(prop,'rate')
                    s.optimizerSettings.uncOptimizerSettings.lineSearchSettings.rate = value;
                elseif strcmp(prop,'maxiter')
                    s.optimizerSettings.maxIter = value;
                elseif strcmp(prop,'printing')
                    s.optimizerSettings.shallPrint = value;
                elseif strcmp(prop,'printMode')
                    s.optimizerSettings.printMode = value;
                elseif strcmp(prop,'monitoring')
                    s.optimizerSettings.monitoringDockerSettings.showOptParams = value;
                elseif strcmp(prop,'monitoring_interval')
                    s.optimizerSettings.monitoringDockerSettings.refreshInterval = value;
                elseif strcmp(prop,'ptype')
                    s.optimizerSettings.monitoringDockerSettings.scale = value;
                elseif strcmp(prop,'plotting')
                    s.optimizerSettings.monitoringDockerSettings.shallDisplayDesignVar = value;
                elseif strcmp(prop,'showBC')
                    s.optimizerSettings.monitoringDockerSettings.shallShowBoundaryConditions = value;
                elseif strcmp(prop,'isDesignVariableFixed')
                    s.designVarSettings.isFixed = value;
                end
            end
            if strcmp(s.designVarSettings.type,'MicroParams')
                s.designVarSettings.creatorSettings.homogSettings = s.homogenizedVarComputerSettings;
            end

            obj.exportFile(s);
        end

    end

    methods (Access = private)

        function init(obj,oldSettings)
            obj.oldSettings = oldSettings;
            obj.propList = fieldnames(oldSettings);
            obj.nProps = length(obj.propList);
        end

        function exportFile(obj,s)
            obj.getFilePath();
            str = jsonencode(s);
            fid = fopen([obj.filePath,'.json'],'w+');
            fprintf(fid,str);
            fclose(fid);
        end

        function getFileName(obj)
            old = obj.oldSettings.case_file;
            new = old;
            %new = [old '.json'];
            obj.fileName = new;
        end

        function getFilePath(obj)
            obj.getFileName();
            obj.filePath = fullfile('.','tests','Applications',obj.fileName);
        end

        function s = translateLevelSetCreator(obj,s)
            prop = obj.oldSettings.levelSetDataBase;
            if ~isempty(prop)
                fields = fieldnames(prop);
                for i = 1:length(fields)
                    field = fields{i};
                    value = prop.(field);
                    s.designVarSettings.creatorSettings.(field) = value;
                end
            end
        end

        function s = translateHomogenizedVariablesComputer(obj,s)
            old = obj.oldSettings;
            switch s.homogenizedVarComputerSettings.type
                case 'ByInterpolation'
                    if isprop(old,'materialInterpolation')
                        s.homogenizedVarComputerSettings.interpolation = old.materialInterpolation;
                    end
                    if isprop(old,'material')
                        s.homogenizedVarComputerSettings.typeOfMaterial = old.material;
                    end
                    if isprop(old,'TOL')
                        s.homogenizedVarComputerSettings.constitutiveProperties = old.TOL;
                    end
                case 'ByVademecum'
                    if isprop(old,'vademecumFileName')
                        s.homogenizedVarComputerSettings.fileName = old.vademecumFileName;
                    end
            end
        end

    end

end