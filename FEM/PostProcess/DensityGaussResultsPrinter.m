classdef DensityGaussResultsPrinter < ResultsPrinter
    
    properties
        fieldName = 'RegularizedDensity';
    end
    
    methods (Access = public)
        
        function obj = DensityGaussResultsPrinter(fileID,fileName,nsteps,gaussDescriptor,etype,ptype,ngaus,ndim,posgp,results,iter)
            obj.init(fileID,fileName,nsteps,gaussDescriptor,etype,ptype,ngaus,ndim,posgp,results,iter)
            obj.print()
        end
    end
    
    
    methods (Access = protected)
        
        function printHeader(obj)
           obj.printGaussPointsHeader()
        end
        
        function printResults(obj)
            dens = obj.results; 
            iS = obj.istep;
            ScalarPrinter(obj.fileID,dens, obj.fieldName,iS,'OnGaussPoints');            
        end
    
    end
    
    methods (Access = private)
        
        function printGaussPointsHeader(obj)
            iD = obj.fileID;
            fprintf(iD,'GaussPoints "%s" Elemtype %s\n',obj.gaussDescriptor,obj.etype);
            fprintf(iD,'Number of Gauss Points: %.0f\n',obj.ngaus);
            fprintf(iD,'Nodes not included\n');
            fprintf(iD,'Natural Coordinates: given\n');
            for igaus = 1:obj.ngaus
                for idime = 1:obj.ndim
                    fprintf(iD,'%12.5d ',obj.posgp(igaus,idime));
                end
                fprintf(iD,'\n');
            end
            fprintf(iD,'End GaussPoints\n');
        end
        
    end
    
    
end