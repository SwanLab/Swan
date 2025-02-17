classdef DensityPrinterForPerimeter < handle
    
    properties (Access = private)
        inputFile
        mesh
        perimeter
        iter
    end
    
    methods (Access = public)
       
        function obj = DensityPrinterForPerimeter(cParams)
           obj.init(cParams)
        end
        
        function print(obj)
            printer = obj.createPrinter();            
            d.cost.shapeFunctions{1} = obj.perimeter;
            d.constraint.shapeFunctions = [];
            quad = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature('CONSTANT');            
            d.quad = quad;
            printer.print(obj.iter,d);               
        end
        
    end
    
    methods (Access = private)
       
        function init(obj,cParams)
            obj.inputFile = cParams.inputFile;
            obj.mesh      = cParams.mesh;
            obj.perimeter = cParams.perimeter;
            obj.iter      = cParams.iter;
        end
        
        function p = createPrinter(obj)
            sP.mesh    = obj.mesh;
            sP.outName = obj.inputFile;
            sP.pdim    = '2D';
            sP.ptype   = 'TRIANGLE';
            p = PostProcessDataBaseCreator(sP);
            s = p.getValue();
            s.cost.shapeFunctions{1} = obj.perimeter;
            s.constraint.shapeFunctions = [];
            type = 'ShapeFunction';
            p = Postprocess(type,s);            
        end                     
        
    end    
    
end