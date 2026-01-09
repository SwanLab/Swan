classdef SubUnfittedMeshSplitter < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        subMesh
    end
       
    methods (Access = public, Static)
        
        function obj = create(cParams)            
            f = SubUnfittedMeshSplitterFactory;
            obj = f.create(cParams);
        end
        
    end
    
    methods (Access = public)
        
        function mI = split(obj)
          
        
    end
    
    methods (Access = protected)
        
        function init(obj,cParams)
            obj.subMesh = cParams.subMesh;
        end    
           
    end
    
    methods (Access = private, Static)
        
            
        
    end
    
    methods (Access = protected, Static)
        createSubMeshComponent(obj)        
    end
    
end