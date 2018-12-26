classdef LevelSetResultsPrinter < ResultsPrinter
    
    properties
        fieldName = 'LevelSet';
    end
    
    
    methods (Access = public)
        
        function obj = LevelSetResultsPrinter(fileID,fileName,nsteps,gaussDescriptor,etype,ptype,ngaus,ndim,posgp,results,iter)
            obj.init(fileID,fileName,nsteps,gaussDescriptor,etype,ptype,ngaus,ndim,posgp,results,iter)
            obj.print()
        end
    end
    
    
    methods (Access = protected)
        
        function printHeader(obj)
        end
        
        function printResults(obj)
            ls = obj.results; 
            iS = obj.istep;
            ScalarPrinter(obj.fileID,ls,obj.fieldName,iS,'OnNodes');            
        end
    
    end
    

end