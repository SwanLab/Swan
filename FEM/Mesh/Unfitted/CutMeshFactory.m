classdef CutMeshFactory < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public, Static)
        
        function c = create(cParams)
            s.coord  = cParams.backgroundMesh.coord;
            s.connec = cParams.backgroundMesh.connec(cParams.cParams,:);
            s.kFace  = cParams.backgroundMesh.kFace;
            backgroundCutMesh = Mesh(s);            

            switch cParams.backgroundMesh.type
                case 'LINE'
                    s.backgroundMesh = backgroundCutMesh;
                    s.cutCells       = cParams.cutCells;
                    s.levelSet       = cParams.levelSet;
                    c = CutMeshProvisionalLine(s);
                case 'TRIANGLE'
                    s.backgroundMesh = backgroundCutMesh;
                    s.cutCells       = cParams.cutCells;
                    s.levelSet       = cParams.levelSet;
                    c = CutMeshComputerProvisional(s);
                                                           
                case 'QUAD'
                    s.backgroundMesh = backgroundCutMesh;
                    s.cutCells       = cParams.cutCells;
                    s.levelSet       = cParams.levelSet;
                    s.lastNode       = max(cParams.backgroundMesh.connec(:));
                    c = CutMeshProvisionalQuadrilater(s);
                otherwise
                    c = CutMeshProvisionalOthers(cParams);                
            end

                
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            
        end
        
    end
    
end