classdef HomogenizedMicrostructureInterpolator < Material
    
    properties (Access = private)
        fileName
        sMesh
        Ctensor
        microParams
    end
    
    methods (Access = public)
        
        function obj = HomogenizedMicrostructureInterpolator(cParams)
            obj.init(cParams);
            [mx,my,C] = obj.loadVademecum();
            obj.createStructuredMesh(mx,my);
            obj.createCtensorFunction(C);
        end

        function C = evaluate(obj,xV)
            C = obj.computeValues(xV);
        end

        function dCm = evaluateDerivative(obj,xV)
           obj.microParams = x;
           dCm{iVar} = obj.createMaterial(x);
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
           obj.fileName    = cParams.fileName;
           obj.microParams = cParams.microParams;
        end

         function [mxV,myV,C] = loadVademecum(obj)
            fName = [obj.fileName,'WithAmplificators'];
            matFile   = [fName,'.mat'];
            file2load = fullfile('Vademecums',matFile);
            v = load(file2load);
            var = v.d;
            mxV = var.domVariables.mxV;
            myV = var.domVariables.myV;
             for imx = 1:length(mxV)
                 for imy = 1:length(myV)
                     C(:,:,imx,imy) = var.variables{imx,imy}.('Ctensor');
                 end
             end
         end 

        function createStructuredMesh(obj,mxV,myV)
            s.x = mxV;
            s.y = myV;
            m = StructuredMesh(s);
            obj.sMesh = m;
        end

        function  createCtensorFunction(obj,C)
            m = obj.sMesh.mesh;
             for i = 1:size(C,1)
                 for j = 1:size(C,2)
                     Cij = squeeze(C(i,j,:,:));
                     CijF = LagrangianFunction.create(m, 1, 'P1');
                     CijF.fValues  = Cij(:);
                     obj.Ctensor{i,j} = CijF;
                 end
             end
        end
        
        function C = computeValues(obj,xV)
            [mL,cells] = obj.obtainLocalCoord(xV);
            nGaus = size(xV,2);
            nElem = obj.microParams{1}.mesh.nelem;
            nStre = size(obj.Ctensor,1); 
          %  nDofs = size(mL,2);
            C  = zeros(nStre,nStre,nGaus,nElem);
            for i = 1:nStre
                for j = 1:nStre
                    Cv = obj.Ctensor{i,j}.sampleFunction(mL,cells);
                    Cij(1,1,:,:) = reshape(Cv,nGaus,[]);
                    C(i,j,:,:)   = Cij(1,1,:,:);
                end
            end
        end

        function [mL,cells] = obtainLocalCoord(obj,xV)
            mx = obj.microParams{1};
            my = obj.microParams{2};
            mxG = mx.evaluate(xV);
            myG = my.evaluate(xV);
            mG(:,1) = mxG(:);
            mG(:,2) = myG(:);
            [mL,cells] = obj.sMesh.obtainLocalFromGlobalCoord(mG);
        end
        
    end
    
    methods (Access = private, Static)
        
        function xR = reshapeDesignVariable(x)
            xV = x.value;
            nV = x.nVariables;
            nx = length(xV)/nV;
            xR = cell(nV,1);
            for ivar = 1:nV
                i0 = nx*(ivar-1) + 1;
                iF = nx*ivar;
                xR{ivar} = xV(i0:iF);
            end
        end
        
    end
    
end