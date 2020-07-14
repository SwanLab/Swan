classdef UnfittedBoundaryMesh < handle
    
    properties (Access = public)
       meshes 
    end
    
    properties (Access = private)
        
    end
    
    properties (Access = private)
       ndim 
       boundaryMesh
    end
    
    methods (Access = public)
        
        function obj = UnfittedBoundaryMesh(cParams)
            obj.init(cParams)            
        end
        
        function compute(obj,ls)
                m = obj.boundaryMesh;
                sides = 2;
                nboxFaces = sides*obj.ndim;
                isBoxFaceMeshActive = false([1 nboxFaces]);
                iFace = 0;
                for idime = 1:obj.ndim
                    for iside = 1:sides
                        iFace = iFace + 1;
                        bMesh = m{iFace};
                        nodesInBoxFace = bMesh.nodesInBoxFaces;
                        s.backgroundMesh = bMesh.mesh;
                        s.isInBoundary   = true;                        
                        cParams = SettingsMeshUnfitted(s);
                        boxFaceMesh = UnfittedMesh(cParams);
                        
                        
                        lsBoxFace = ls(nodesInBoxFace);
                        if any(sign(lsBoxFace)<0)
                            boxFaceMesh.compute(lsBoxFace);
                            isBoxFaceMeshActive(iFace) = true;
                        end                        
                        
                        
                        u.boxFaceMeshes       = boxFaceMesh;
                        u.nodesInBoxFaces     = nodesInBoxFace;
                        u.isBoxFaceMeshActive = isBoxFaceMeshActive;

                        obj.meshes{iFace} = u;
                    end
                end           
        end
           
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
          obj.ndim         = cParams.ndim;
          obj.boundaryMesh = cParams.boundaryMesh;                      
        end
        
    end
    
end