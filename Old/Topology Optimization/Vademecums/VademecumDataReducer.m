classdef VademecumDataReducer < handle
    
    properties (Access = private)
       fileName
       fileNameReduced
       fullVademecum
       reducedVademecum
    end
    
    methods (Access = public)
        
        function obj = VademecumDataReducer()
            obj.init();
            obj.loadFullVademecum();
            obj.createReducedVademecum();           
            obj.saveReducedVademecum();
        end
        
    end
        
    methods (Access = private)
        
        function init(obj)
            obj.fileName = 'SuperEllipseQ2';%'SuperEllipseQOptAnalytic';%'SuperEllipseQOptAnalytic';%'SuperEllipseQMax';%'SuperEllipseQ2';            
            obj.fileNameReduced = [obj.fileName,'Reduced'];
        end
        
        function loadFullVademecum(obj)
            matFile   = [obj.fileName,'.mat'];
            file2load = fullfile('/media/alex/MyPassport/Vademecum/',matFile);
            v = load(file2load);            
            obj.fullVademecum = v.d;
        end
        
        function createReducedVademecum(obj)
            fv = obj.fullVademecum;
            rv.fileName = obj.fileNameReduced;
            rv.domVariables = fv.domVariables;
            rv.outPutPath   = fv.outPutPath;
            nx = length(fv.domVariables.mxV);
            ny = length(fv.domVariables.mxV);
            for im = 1:nx
                for jm = 1:ny
                    fvIJ = fv.variables{im,jm};
                    rv.variables{im,jm}.volume  = fvIJ.volume;
                    rv.variables{im,jm}.Ctensor = fvIJ.Ctensor;
                    rv.variables{im,jm}.integrationVar = fvIJ.integrationVar;
                end
            end  
            obj.reducedVademecum = rv;
        end
        
        function saveReducedVademecum(obj)
            d = obj.reducedVademecum;
            matFile   = [obj.fileNameReduced,'.mat'];
            file2save = fullfile('Topology Optimization','Vademecums',matFile);            
            save(file2save,'d');
        end
         
        
    end
        
end
