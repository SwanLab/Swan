classdef Printer < handle
    
    properties (Access = private)
       mesh 
    end
    
    
    properties (Access = protected)
       quadrature        
       PostProcess
    end
    
    methods (Access = public)
        
        function print(obj,variable,outname)
            res = obj.createResultsInputStructure(variable,outname);
            m = obj.mesh;
            obj.PostProcess.print(m,res)
        end
        
    end
    
    methods (Access = protected)
        
        function init(obj,quad,mesh)
            obj.quadrature = quad;
            obj.mesh = mesh;
            obj.createPostProcess()
        end
        
    end
    
    methods (Abstract,Access =  protected, Static)
     createPostProcess(obj)        
     createResultsInputStructure(obj,dens,outname)
    end
    
end
