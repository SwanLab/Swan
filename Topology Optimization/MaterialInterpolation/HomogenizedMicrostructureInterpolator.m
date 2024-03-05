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

        function C = obtainTensor(obj)
          s.operation = @(xV) obj.evaluate(xV);
          C = DomainFunction(s);            
        end
        
        function dC = obtainTensorDerivative(obj)
          s.operation = @(xV) obj.evaluateGradientM1(xV);
          s.ndimf = 1;
          dC{1} =  DomainFunction(s);
          s.operation = @(xV) obj.evaluateGradientM2(xV);
          s.ndimf = 1;
          dC{2} =  DomainFunction(s);
        end             

        function setDesignVariable(obj,x)
            obj.microParams = x;
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
        
        function C = evaluate(obj,xV)
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

        function dC = evaluateGradientM1(obj,xV)
            dC = obj.evaluateGradient(xV,1);
        end

        function dC = evaluateGradientM2(obj,xV)
            dC = obj.evaluateGradient(xV,2);
        end        

        function dCt = evaluateGradient(obj,xV,dir)
            [mL,cells] = obj.obtainLocalCoord(xV);
            nGaus = size(xV,2);
            nElem = obj.microParams{1}.mesh.nelem;
            nDim  = obj.microParams{1}.mesh.ndim;
            nStre = size(obj.Ctensor,1);
            %  nDofs = size(mL,2);
            dC  = zeros(nDim,nStre,nStre,nGaus,nElem);
            for i = 1:nStre
                for j = 1:nStre
                    dCv = obj.Ctensor{i,j}.sampleGradient(mL,cells);
                    dCv = squeezeParticular(dCv,1);
                    dCij(:,1,1,:,:)  = reshape(dCv,nDim,nGaus,nElem);
                    dC(:,i,j,:,:) = dCij;
                end
            end     
            dCt = squeezeParticular(dC(dir,:,:,:,:),1);
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