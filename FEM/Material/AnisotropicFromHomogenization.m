classdef AnisotropicFromHomogenization < Material
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        vadVariables
        interpolator        
        Ctensor 
        sMesh
    end
    
    properties (Access = private)
        microParams
        fileName
    end
    
    methods (Access = public)
        
        function obj = AnisotropicFromHomogenization(cParams)
            obj.init(cParams)
            obj.loadVademecum();    
            obj.createStructuredMesh();
            obj.obtainValues();      
            obj.createInterpolator();
        end
        
        function C = evaluate(obj,xV)
            C = obj.computeValues(xV);
            
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
           obj.microParams = cParams.microParams;
           obj.fileName = cParams.fileName;           
        end
        
        function loadVademecum(obj)
            fName = [obj.fileName,'WithAmplificators'];
            matFile   = [fName,'.mat'];
            file2load = fullfile('Vademecums',matFile);
            v = load(file2load);
            obj.vadVariables = v.d;          
        end
        
        function obtainValues(obj)
            var = obj.vadVariables.variables;
            mxV = obj.vadVariables.domVariables.mxV;
            myV = obj.vadVariables.domVariables.mxV;
            for imx = 1:length(mxV)
                for imy = 1:length(myV)
                    C(:,:,imx,imy) = var{imx,imy}.('Ctensor');
                end
            end   
            
            for i = 1:size(C,1)
                for j = 1:size(C,2)
                    Cij = squeeze(C(i,j,:,:));
                    s.fValues = Cij(:);
                    s.mesh    = obj.sMesh.mesh;
                    s.ndim    = 1;
                    CijF      = P1Function(s);
                    obj.Ctensor{i,j} = CijF;
                end
            end
            
        end         

        function createStructuredMesh(obj)
            s.x = obj.vadVariables.domVariables.mxV;
            s.y = obj.vadVariables.domVariables.myV;
            m = StructuredMesh(s); 
            obj.sMesh = m;
        end
        
        function createInterpolator(obj)
            s.mesh = obj.sMesh.mesh;
            obj.interpolator = Interpolator(s);            
        end

        function [mxG,myG] = computeGlobalEvaluationPoints(obj,xV)
            mx = obj.microParams{1};
            my = obj.microParams{2};
            mxG = mx.evaluate(xV);
            myG = my.evaluate(xV);
        end

        function [xL,cells] = computeLocalEvaluationPoints(obj,xG1,xG2)
            


            s.mesh     = obj.sMesh;
            s.points.x = (xG1(:));
            s.points.y = (xG2(:));
            cellFinder = CellFinderInStructuredMesh(s);
            xL    = cellFinder.naturalCoord;
            cells = cellFinder.cells;
%             nC = size(xG1,3);            
%             xG1 = permute(xG1,[2 1 3]);
%             xG1 = reshape(x,[],nC);
%             nodes = obj.cellFinder.cells;
%             [nelem,nnode] = size(nodes);
%             xNode = x(nodes(:),:);
%             xN = reshape(xNode,nelem,nnode,nC);            
        end        
        
        
     function C = computeValues(obj,xV)
            [mxG,myG] = obj.computeGlobalEvaluationPoints(xV);
            [mL,cells] = obj.computeLocalEvaluationPoints(mxG,myG);                       
         %   mL = reshape(mL,[2,3,9604]);
            nStre = size(obj.Ctensor,1);            
           % C  = zeros(nStre,nStre,mx.nDofs);
            for i = 1:nStre
                for j = 1:nStre
                    Cij = obj.Ctensor{i,j}.evaluate(mL);  
                    C(i,j,:) = Cij;
                end
            end

            %nVar = length(obj.microParams);            
            %dC = zeros(nStre,nStre,mx.nDofs,nVar);

            % for i = 1:nStre
            %     for j = 1:nStre
            %         for iVar = 1:nVar
            %             dC(i,j,:,iVar) = dc(:,iVar);
            %         end
            %     end
            % end  
            %dCt{1} = dC(:,:,:,1);
            %dCt{2} = dC(:,:,:,2);
        end                  
        
    end
    
end