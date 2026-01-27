classdef LinearInterpolator < handle

    properties (Access = public)
        f
        df
    end

    properties (Access = private)
        mx,my,Ch
        sMesh
    end

    properties (Access = private)

    end

    methods (Access = public)

        function obj = LinearInterpolator()
             d = load('PaperDehomog/MicroStructureOptimization/OfflineRectangularChomog.mat');
            obj.mx = d.mx;
            obj.my = d.my;
            obj.Ch = d.Ch;
            obj.createStructuredMesh(obj.mx,obj.my);
            obj.createCtensorFunction(C);
        end

    end

    methods (Access = private)

        function init(obj)

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
                     for k = 1:size(C,3)
                         for l = 1:size(C,4)
                     Cij = squeeze(C(i,j,:,:));
                     CijF = LagrangianFunction.create(m, 1, 'P1');
                     CijF.fValues  = Cij(:);
                     obj.f{i,j} = CijF;
                         end
                     end
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
                    Cv = obj.f{i,j}.sampleFunction(mL,cells);
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

end