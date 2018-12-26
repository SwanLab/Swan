classdef DensityResultsPrinter < ResultsPrinter
    
    properties
        fieldName = 'Density';
    end
    
    methods (Access = public)
        
        function obj = DensityResultsPrinter(fileID,fileName,nsteps,gaussDescriptor,etype,ptype,ngaus,ndim,posgp,results,iter)
            obj.init(fileID,fileName,nsteps,gaussDescriptor,etype,ptype,ngaus,ndim,posgp,results,iter)
            obj.print()
        end
    end
    
    
    methods (Access = protected)
        
        function printHeader(obj)
        end
        
        function printResults(obj)
            dens = obj.results; 
            iS = obj.istep;
            ScalarPrinter(obj.fileID,dens, obj.fieldName,iS,'OnNodes');            
        end
    
    end
    
    
end