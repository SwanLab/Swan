classdef InnerCutMeshSplitter < SubUnfittedMeshSplitter
    
    
    methods (Access = public)
        
        function obj = InnerCutMeshSplitter(cParams)
            obj.init(cParams)
        end
        
    end
    
    methods (Access = protected, Static)
        
        function sM = createSubMeshComponent(m,nodes)

        end        
        
    end
    
end