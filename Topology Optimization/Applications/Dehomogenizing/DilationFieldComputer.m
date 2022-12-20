classdef DilationFieldComputer < handle
    
    properties (Access = private)
       dilation 
    end
    
    properties (Access = private)
       theta 
       mesh
    end
    
    methods (Access = public)
        
        function obj = DilationFieldComputer(cParams)
            obj.init(cParams)
        end
        
        function d = compute(obj)
            obj.computeDilationField();
            d = obj.dilation; 
        end
        
        function plot(obj)
            figure()
            s.mesh  = obj.mesh;
            s.field = obj.dilation;
            n = NodalFieldPlotter(s);
            n.plot();
            shading interp            
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.theta = cParams.theta;
            obj.mesh  = cParams.mesh;
        end
               
        function computeDilationField(obj)
            s.fGauss = obj.computeThetaGradient();
            s.mesh   = obj.mesh;
            varProb  = MinimumGradFieldWithVectorInL2(s);
            r = varProb.solve();
            obj.dilation = r;
        end
        
        function gradT = computeThetaGradient(obj)
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature('CUBIC');
            m.type = obj.mesh.type;
            int = Interpolation.create(m,'LINEAR');
            int.computeShapeDeriv(q.posgp);
            s.mesh = obj.mesh;
            g = Geometry.create(s);
            g.computeGeometry(q,int);
            grad = g.dNdx;
            nodes = obj.mesh.connec;
            f = obj.theta;
            gradF = zeros(obj.mesh.ndim,q.ngaus,obj.mesh.nelem);
            for igaus = 1:q.ngaus
                for kNode = 1:obj.mesh.nnodeElem
                    nodeK = nodes(:,kNode);
                    fI = f(nodeK);
                    for idim = 1:obj.mesh.ndim
                        dN = squeeze(grad(idim,kNode,:,igaus));
                        gF = squeeze(gradF(idim,igaus,:));
                        gradF(idim,igaus,:) = gF + fI.*dN;%bsxfun(@times,dN,fI);
                    end
                end
            end
            gradT = zeros(size(gradF));
            gradT(1,:,:) = -gradF(2,:,:);
            gradT(2,:,:) = gradF(1,:,:);

            alpha = obj.theta;
            a1(:,1) = cos(alpha);
            a1(:,2) = sin(alpha);

            a2(:,1) = -sin(alpha);
            a2(:,2) = cos(alpha);

            dN = g.dNdx;
            N  = int.shape;

            Q = [0 1;-1 0];
            nNode = obj.mesh.nnodeElem;
            nDim = obj.mesh.ndim;
            gradT2 = zeros(size(gradF));
            for igaus = 1:q.ngaus
                for jNode = 1:nNode
                     nodeJ = nodes(:,jNode);                                            
                    for kNode = 1:nNode
                        nodeK = nodes(:,kNode);                        
                        for rDim = 1:nDim
                            for lDim = 1:nDim
                                for mDim = 1:nDim
                                    dNkl = squeeze(dN(lDim,kNode,:,igaus));
                                    Nj   = N(jNode,igaus);
                                    Qlr = Q(lDim,rDim);
                                    a1kr = a1(nodeK,rDim);
                                    a1jm = a1(nodeJ,mDim);
                                    a2kr = a2(nodeK,rDim);
                                    a2jm = a2(nodeJ,mDim);
                                    gV(1,1,:) = dNkl.*Nj.*Qlr.*(a1kr.*a2jm - a2kr.*a1jm);
                                    gradT2(mDim,igaus,:) = gradT2(mDim,igaus,:) + gV;
                                end
                            end
                        end
                    end
                end
            end


            %                    nodeK = nodes(:,kNode);
            %                    fI = f(nodeK);
            %                    for idim = 1:obj.mesh.ndim
            %                     dN = squeeze(grad(idim,kNode,:,igaus));
            %                     gF = squeeze(gradF(idim,igaus,:));
            %                     gradF(idim,igaus,:) = gF + fI.*dN;%bsxfun(@times,dN,fI);
            %                    end
            %                 end
            gradT = gradT2;
        end




    end

    
end