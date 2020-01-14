classdef DensityPrinterForPerimeter < handle
    
    properties (Access = private)
        inputFile
        mesh
        density
        iter
    end
    
    methods (Access = public)
       
        function obj = DensityPrinterForPerimeter(cParams)
           obj.init(cParams)
        end
        
        function print(obj)
            printer = obj.createPrinter();
            d.x  = obj.density;
            printer.print(obj.iter,d);               
        end
        
    end
    
    methods (Access = private)
       
        function init(obj,cParams)
            obj.inputFile = cParams.inputFile;
            obj.mesh      = cParams.mesh;
            obj.density   = cParams.density;
            obj.iter      = cParams.iter;
        end
        
        function p = createPrinter(obj)
            sP.mesh    = obj.mesh;
            sP.outName = obj.inputFile;
            sP.pdim    = '2D';
            sP.ptype   = 'TRIANGLE';
            p = PostProcessDataBaseCreator(sP);
            s = p.getValue();
            type = 'Density';
            p = Postprocess(type,s);            
        end                     
        
    end    
    
end