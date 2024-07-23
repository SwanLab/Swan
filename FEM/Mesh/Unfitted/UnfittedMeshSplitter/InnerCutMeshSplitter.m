classdef InnerCutMeshSplitter < SubUnfittedMeshSplitter
    
    
    methods (Access = public)
        
        function obj = InnerCutMeshSplitter(cParams)
            obj.init(cParams)
        end
        
    end
    
    methods (Access = protected, Static)
        
        function sM = createSubMeshComponent(m,nodes)
            s.connec  = m.mesh.connec(nodes,:);
            s.coord   = m.mesh.coord;
            s.mesh                  = Mesh.create(s);
            s.xCoordsIso            = m.xCoordsIso(:,:,nodes);
            s.cellContainingSubcell = m.cellContainingSubcell(nodes,:);
            sM = InnerCutMesh(s);
        end        
        
    end
    
end