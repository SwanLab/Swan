classdef ImageTakerFor3DShape < handle
    
    properties (Access = private)
        uMesh
        outDir
    end
    
    methods (Access = public)
        
        function obj = ImageTakerFor3DShape(cParams)
            obj.init(cParams)            
        end
        
        function takeImages(obj)
            obj.plotUnfittedMesh();
            obj.takeIsoFrontImage();
            obj.takeBackImage();
            obj.takeFrontImage();
            obj.takeIsoBackImage();
        end        
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.uMesh  = cParams.unfittedMesh;
            obj.outDir = cParams.outDir;
        end
        
        function takeIsoFrontImage(obj)
            fName = 'IsoFrontImage';
            az = 120;
            el = 30;
            fName = [obj.outDir,fName];
            view(az,el)
            obj.takeImage(fName)
        end
        
        function takeIsoBackImage(obj)
            fName = 'IsoBackImage';
            az = 240;
            el = 30;
            fName = [obj.outDir,fName];
            view(az,el)
            obj.takeImage(fName)
        end
        
        function takeFrontImage(obj)
            fName = 'YZfrontImage';
            az = 90;
            el = 0;
            fName = [obj.outDir,fName];
            view(az,el)
            obj.takeImage(fName)
        end
        
        function takeBackImage(obj)
            fName = 'YZbackImage';
            az = 90;
            el = 180;
            fName = [obj.outDir,fName];
            view(az,el)
            obj.takeImage(fName)
        end    
        
        function plotUnfittedMesh(obj)
            f = figure();
            f.WindowState = 'maximized';
            obj.uMesh.plot();            
        end
        
        
    end
        
    methods (Access = private, Static)
                
        function takeImage(name)
            print(name,'-dpng')
        end
        
    end
    
end