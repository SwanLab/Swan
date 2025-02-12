classdef OptimizationMetricsPrinter < handle
    
    properties (Access = protected)
        optimizer
    end
    
    properties (Access = private)
        filePath
        folderPath
        cost
        constraint
    end
    
    methods (Access = protected, Abstract)
        printConvergenceVariables(obj,fid)
    end
    
    methods (Access = public)
        
        function obj = OptimizationMetricsPrinter(cParams)
            obj.init(cParams);
        end
        
        function print(obj,nIter,iStep)
            fid = fopen(obj.filePath,'at');
            
            fprintf(fid,'-----------------------------------------------------------------------------------------------\n');
            fprintf(fid,'\n');
            fprintf(fid,'Iteration: %i \n',nIter);
            fprintf(fid,'Nstep: %i \n',iStep);
            fprintf(fid,'cost %f \n',obj.cost.value);
            
            for i = 1:length(obj.cost.shapeFunctions)
                fprintf(fid,strcat('-cost ',num2str(i),': %f \n'),obj.cost.shapeFunctions{i}.value);
            end
            
            fprintf(fid,'constraint: %f \n',obj.constraint.value);
            for i = 1:length(obj.constraint.shapeFunctions)
                fprintf(fid,strcat('-constraint ',num2str(i),': %f \n'),obj.constraint.shapeFunctions{i}.value);
            end
            
           % obj.printConvergenceVariables(fid);
            
            fprintf(fid,'\n');
            
            fclose(fid);
        end
        
        function printFinal(obj)
           fid = fopen(obj.filePath,'at');
           fprintf(fid,'\n');           
           fprintf(fid,'-----------------------------------------------------------------------------------------------\n');
           fprintf(fid,'\n');
           fprintf(fid,'cost %f \n',obj.cost.computeNonNormalizedValue());          
           fclose(fid);
        end
        
    end
    
    methods (Access = protected)
        
        function init(obj,cParams)
            obj.createPaths(cParams);
            obj.optimizer = cParams.optimizer;
            obj.cost = cParams.cost;
            obj.constraint = cParams.constraint;
            obj.createFolderIfNeeded();
            obj.createFile();
        end
        
    end
    
    methods (Access = private)
        
        function createPaths(obj,cParams)
            fileName = cParams.fileName;
            obj.folderPath = fullfile('Output',fileName);
            obj.filePath = fullfile(obj.folderPath,strcat(fileName,'.txt'));
        end
        
        function createFile(obj)
            fid = fopen(obj.filePath,'wt');
            fclose(fid);
        end
        
        function createFolderIfNeeded(obj)
            if ~exist(obj.folderPath,'dir')
                mkdir(obj.folderPath);
            end
        end
        
    end
    
end