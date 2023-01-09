classdef P1DiscontinuousFunction < FeFunction
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        connec
        fNodes
        interpolation
        type
    end
    
    properties (Access = private)
    end
    
    methods (Access = public)
        
        function obj = P1DiscontinuousFunction(cParams)
            obj.init(cParams)
            obj.createInterpolation();
        end

        function fxV = evaluate(obj, xV)
            % Goal: develop this function
            func = obj.fValues;
            obj.interpolation.computeShapeDeriv(xV);
            shapes = obj.interpolation.shape;
            nNode  = size(shapes,1);
            nGaus  = size(shapes,2);
            nF     = size(func,1);
            nElem  = size(func,3);
            fxV = zeros(nF,nGaus,nElem);
            for kNode = 1:nNode
                shapeKJ = shapes(kNode,:,:);
                fKJ     = func(:,kNode,:);
                f = bsxfun(@times,shapeKJ,fKJ);
                fxV = fxV + f;
            end
        end
        
        function fFine = refine(obj,m,mFine)
         %   mFineD = mFine.createDiscontinuousMesh();
            f = squeeze(obj.fValues); 
            f = f(:);
            fEdges = obj.computeFunctionInEdges(m,f);                   
            fAll  = [f;fEdges]; 
           
            
            s.type    = mFine.type;
            s.connec  = mFine.connec;
            s.fValues = fAll;
            fP1 = P1Function(s);
            
            s.mesh   = mFine;
            s.connec = mFine.connec;
            s.type   = mFine.type;
            sP.origin = 'P1';
            sP.x = fP1;

            p = ProjectorToP1discont(s);
            fFine = p.project(sP);            
%             
%             for iNode = 1:m.nnodeElem  
%                 fV(1,iNode,:) = fAll(mFine.connec(:,iNode));
%             end
%           %  ndims = size(fAll, 2);
%           %  nNodesDisc = m.nnodeElem*m.nelem;
%           %  coordD = reshape(obj.xFE.fValues, [ndims, nNodesDisc])';
%             s.type = obj.type;
%             s.connec = obj.connec;
%             s.fValues = obj.coord;
%             coordP1 = P1Function(s);
% 
%             s.type    = mFine.type;
%             s.connec  = mFineD.connec;
%             s.fValues = fV;     
%             fFine = P1DiscontinuousFunction(s);
        end

        function plot(obj, m)
            ndims   = size(obj.fValues, 1);
            nelem   = size(obj.connec, 1);
            nnodeEl = size(obj.connec, 2);
            mD = m.createDiscontinuousMesh();
            fD = reshape(obj.fValues, [ndims, nelem*nnodeEl])';
            x = mD.coord(:,1);
            y = mD.coord(:,2);
            figure()
            for idim = 1:obj.ndimf
                subplot(1,obj.ndimf,idim);
                z = fD(:,idim);
                a = trisurf(mD.connec,x,y,z);
                view(0,90)
    %             colorbar
                shading interp
                a.EdgeColor = [0 0 0];
                title(['dim = ', num2str(idim)]);
            end
        end

    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.fValues = cParams.fValues;
            obj.connec  = cParams.connec;
            obj.type    = cParams.type;
            obj.ndimf   = size(cParams.fValues,1);
        end

        function createInterpolation(obj)
            m.type = obj.type;
            obj.interpolation = Interpolation.create(m,'LINEAR');
        end

        function f = computeFunctionInEdges(obj,m,fNodes)
            s.mesh   = m;
            s.fNodes = fNodes;
            eF       = EdgeFunctionInterpolator(s);
            f = eF.compute();
        end        
        
    end
    
end