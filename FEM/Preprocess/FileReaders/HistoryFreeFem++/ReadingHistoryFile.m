classdef ReadingHistoryFile < handle
    
    properties (SetAccess = private, GetAccess = public)
        cost
        lambda
    end
    
    properties (Access = private)
        filePath
        fid
        nIter
    end
    
    methods (Access = public)
        
        function obj = ReadingHistoryFile(d)
            obj.init(d);
            obj.obtainNumberOfIterations();
            obj.obtainHistoryData();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,d)
            obj.filePath = d.filePath;
        end
        
        function obtainNumberOfIterations(obj)
            obj.openFile();
            obj.nIter = 0;
            tline = fgetl(obj.fid);
            while ischar(tline)
                tline = fgetl(obj.fid);
                obj.nIter = obj.nIter + 1;
            end
            obj.nIter = obj.nIter;
            obj.closeFile();
        end
        
        function obtainHistoryData(obj)
            obj.openFile();
            obj.readCost();
            obj.closeFile();
        end
        
        function openFile(obj)
            obj.fid = fopen(obj.filePath);
        end
        
        function closeFile(obj)
            fclose(obj.fid);
        end
        
        
        function readCost(obj)
            for ielem=1:obj.nIter
                line =  obj.readLine();
                obj.cost(ielem) = line(2);
                obj.lambda(ielem) = line(3);
            end
        end
        
        function var = readLine(obj)
            line = fgetl(obj.fid);
            var = str2num(line);
        end
        
        
    end
    
    
    
end