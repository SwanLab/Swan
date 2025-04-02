classdef VademecumAmplificatorAdder < handle
    
    properties (Access = private)
       fileName
       fullFileName
       fileNameReduced
       fullFileNameReduced
       reducedVademecum
       fullVademecum
       fullVademecumWithAmplificators
       reducedVademecumWithAmplificators    
       order
    end
    
    methods (Access = public)
        
        function obj = VademecumAmplificatorAdder()
            obj.init();
            obj.loadReducedVademecum();
            obj.computeFullVademecumWithAmplificators();
            obj.createReducedVademecumWithAmplificators();
            obj.saveReducedVademecum();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.fileName        = 'SuperEllipseQOptAnalytic';%'SuperEllipseQMax';'SuperEllipseQOptAnalytic';'SuperEllipseQ2';
            obj.fileNameReduced = [obj.fileName,'Reduced'];
            obj.fullFileName        = obj.obtainFullFileName(obj.fileName);  
            obj.fullFileNameReduced = obj.obtainFullFileName(obj.fileNameReduced);  
            obj.order = 2:2:8;%[2,4,6,8,10,12,14,16];            
        end
        
        function loadReducedVademecum(obj)
            v = load(obj.fileNameReduced);
            obj.reducedVademecum = v.d;             
        end
        
        function computeFullVademecumWithAmplificators(obj)
            for iorder = 1:length(obj.order)
                iorder
                cParams.fileName = obj.fullFileName;
                cParams.pNorm    = obj.order(iorder);
                vc = VademecumPtensorComputer(cParams);
                vc.compute();
                vd{iorder} = vc.vademecumData;
            end           
            obj.fullVademecumWithAmplificators = vd;
        end        
        
        function createReducedVademecumWithAmplificators(obj)
            obj.obtainAmplificators();
            obj.obtainMonomials();
        end
        
        function obtainAmplificators(obj)
            rv = obj.reducedVademecum;
            nx = length(rv.domVariables.mxV);
            ny = length(rv.domVariables.mxV);
            for im = 1:nx
                for jm = 1:ny
                    rv.variables{im,jm}.Ptensor  = obj.obtainPtensor(im,jm);
                end
            end            
            obj.reducedVademecumWithAmplificators = rv;                        
        end

        function Ptensor = obtainPtensor(obj,im,jm)
            v = obj.fullVademecumWithAmplificators;
            nOrder = length(obj.order);
            Ptensor = cell(nOrder,1);
            for iorder = 1:nOrder
                vi = v{iorder};
                Ptensor{iorder} = vi.variables{im,jm}.Ptensor;
            end            
        end
        
        function obtainMonomials(obj)
            v = obj.fullVademecumWithAmplificators;
            nOrder = length(obj.order);
            monomials = cell(nOrder,1);            
            for iorder = 1:nOrder
                vi = v{iorder};
                monomials{iorder} = vi.monomials;
            end      
            obj.reducedVademecumWithAmplificators.monomials = monomials;                                    
        end
        
        function saveReducedVademecum(obj)
            d = obj.reducedVademecumWithAmplificators;
            matFile   = [obj.fileName,'WithAmplificators','.mat'];
            file2save = fullfile('Topology Optimization','Vademecums',matFile);
            save(file2save,'d');
        end        
        
    end
    
    methods (Access = private, Static)
        
        function f = obtainFullFileName(fileName)
            matFile   = [fileName,'.mat'];
            f = ['/media/alex/MyPassport/Vademecum/',matFile];                        
        end        
        
    end
    
end
