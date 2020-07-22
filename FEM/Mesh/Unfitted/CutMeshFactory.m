classdef CutMeshFactory < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public, Static)
        
        function c = create(cParams)
            bMesh = cParams.backgroundMesh;
            s.coord  = bMesh.coord;
            s.connec = bMesh.connec(cParams.cutCells,:);
            s.kFace  = bMesh.kFace;
            backgroundCutMesh = Mesh(s);            

            switch cParams.backgroundMesh.type
                case 'LINE'
                    s.backgroundMesh    = bMesh;
                    s.backgroundCutMesh = backgroundCutMesh;
                    s.cutCells          = cParams.cutCells;
                    s.levelSet          = cParams.levelSet;
                    c = CutMeshProvisionalLine(s);
                case 'TRIANGLE'
                    s.backgroundMesh    = bMesh;
                    s.backgroundCutMesh = backgroundCutMesh;                    
                    s.cutCells       = cParams.cutCells;
                    s.levelSet       = cParams.levelSet;
                    c = CutMeshComputerProvisional(s);
                                                           
                case 'QUAD'
                    s.backgroundMesh    = bMesh;                    
                    s.backgroundCutMesh = backgroundCutMesh;
                    s.cutCells          = cParams.cutCells;
                    s.levelSet          = cParams.levelSet;
                    s.lastNode          = max(cParams.backgroundMesh.connec(:));
                    c = CutMeshProvisionalQuadrilater(s);
                otherwise
                    s.backgroundMesh    = bMesh;                    
                    s.backgroundCutMesh = backgroundCutMesh;
                    s.cutCells          = cParams.cutCells;
                    s.levelSet          = cParams.levelSet;
                    s.type              = cParams.type;
                    c = CutMeshProvisionalOthers(s);                
            end

                
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            
        end
        
    end
    
end