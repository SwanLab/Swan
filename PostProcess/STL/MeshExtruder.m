classdef MeshExtruder < handle
    properties (Access = private)
        height
        filename
        unfittedMesh
    end
    
    methods (Access = public)
        
        function obj = MeshExtruder(cParams)
            obj.init(cParams)
        end

        function m =  extrude(obj)
            s2g = SwanGiDInterface();
            s2g.extrudeMesh(obj.unfittedMesh, obj.height);
            m = obj.readMsh();
            delete MeshExtruder_File.flavia.msh
            delete MeshExtruder_File.flavia.res
            delete PostProcess/STL/sampleMesh.msh
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.height       = cParams.height;
            obj.unfittedMesh = cParams.unfittedMesh;
            obj.filename     = 'MeshExtruder_File';
        end

        function m = readMsh(obj)
            s.filePath = obj.getOutputFileName();
            mR = MshReader(s);
            m = mR.read();
        end

        function f = getOutputFileName(obj)
            f = [pwd,'/PostProcess/STL/sampleMesh.msh'];
        end

        
    end
    
end