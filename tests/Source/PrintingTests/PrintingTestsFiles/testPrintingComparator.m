classdef testPrintingComparator < handle
    
    properties (Access = protected, Abstract)
        fileOutputName
        iter
    end
    
    methods (Access = protected)
       
        function hasPassed = hasPassed(obj)
            hasMshChanged = obj.compareFile('.msh');
            hasResChanged = obj.compareFile('.res');   
            filesHaveChanged = hasMshChanged || hasResChanged;            
            hasPassed = ~filesHaveChanged;
        end
        
        function selectComputedVar(obj)
        end
        
    end
    
    methods (Access = private)
        
       function sF = obtainSavedPrintedFile(obj,ext)
           out   = obj.fileOutputName;            
           fullOutSavedName = [out,'.flavia',ext];            
           sF = fullfile('tests','PrintingTests',out,fullOutSavedName);
        end
        
        function oF = obtainOutPutFile(obj,ext)
           out   = obj.fileOutputName;            
           fullOutName = [out,num2str(obj.iter),'.flavia',ext];
           oF = fullfile('Output',out,fullOutName);              
        end
        
        function hasChanged = compareFile(obj,ext)
            sF = obj.obtainSavedPrintedFile(ext);
            oF = obj.obtainOutPutFile(ext);    
            hasChanged = FileComparator().areFilesDifferent(sF,oF);
        end        
        
    end

    
end    