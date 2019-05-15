classdef SettingsTranslator < handle
    
    properties (GetAccess = public, SetAccess = private)
        fileName
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
%                 if ~isempty(value)
                    if strcmp(prop,'filename')
                        s.problemData.problemFileName = value;
                    elseif strcmp(prop,'ptype')
                        s.problemData.scale = value;
                    elseif strcmp(prop,'materialInterpolation')
                        s.homogenizedVarComputerSettings.interpolation = value;
                    elseif strcmp(prop,'material')
                        s.homogenizedVarComputerSettings.typeOfMaterial = value;
                    elseif strcmp(prop,'homegenizedVariablesComputer')
                        s.homogenizedVarComputerSettings.type = value;
                    elseif strcmp(prop,'TOL')
                        s.homogenizedVarComputerSettings.constitutiveProperties = value;
                    elseif strcmp(prop,'designVariable')
                        s.designVarSettings.type = value;
                    elseif strcmp(prop,'initial_case')
                        s.designVarSettings.initialCase = value;
                    elseif strcmp(prop,'nsteps')
                        s.incrementalSchemeSettings.nSteps = value;
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
                    elseif strcmp(prop,'cost')
                        for k = 1:length(value)
                            s.costSettings.shapeFuncSettings{k}.type = value{k};
                        end
                    elseif strcmp(prop,'weights')
                        s.costSettings.weights = value;
                    elseif strcmp(prop,'filter')
                        s.costSettings.filterType = value;
                        s.constraintSettings.filterType = value;
                    elseif strcmp(prop,'constraint')
                        for k = 1:length(value)
                            s.constraintSettings.shapeFuncSettings{k}.type = value{k};
                        end
                    elseif strcmp(prop,'optimizer')
                        s.optimizerSettings.type = value;
                    elseif strcmp(prop,'optimizerUnconstrained')
                        s.optimizerSettings.uncOptimizerSettings.type = value;
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
                    elseif strcmp(prop,'plotting')
                        s.optimizerSettings.monitoringDockerSettings.shallDisplayDesignVar = value;
                    elseif strcmp(prop,'showBC')
                        s.optimizerSettings.monitoringDockerSettings.shallShowBoundaryConditions = value;
                    end
%                 end
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
            obj.getFileName();
            json = jsonencode(s);
            fid = fopen(obj.fileName,'w+');
            fprintf(fid,json);
            fclose(fid);
        end
        
        function getFileName(obj)
            old = obj.oldSettings.case_file;
            new = [old '.json'];
            obj.fileName = new;
        end
        
    end
    
end