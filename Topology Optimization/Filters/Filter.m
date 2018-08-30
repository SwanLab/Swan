classdef Filter < handle
    properties
        diffReacProb
        M0
        coordinates
        connectivities
        nnode
        nelem
        npnod
        ngaus
        shape
        x
        x_reg
        P_operator
    end
    
    methods
        function obj = Filter(problemID,scale)
            switch scale
                case 'MACRO'
                    obj.diffReacProb = DiffReact_Problem(problemID);
                case 'MICRO'
                    obj.diffReacProb = DiffReact_Problem_Micro(problemID);
            end
        end        
        function preProcess(obj)
            obj.diffReacProb.preProcess;
            quadrature = Quadrature.set(obj.diffReacProb.geometry.type);
            quadrature.computeQuadrature('LINEAR');
            obj.diffReacProb.element.interpolation_u.computeShapeDeriv(quadrature.posgp)
            obj.diffReacProb.geometry.computeGeometry(quadrature,obj.diffReacProb.element.interpolation_u);
            
            for igauss = 1:quadrature.ngaus
                obj.M0{igauss} = sparse(1:obj.diffReacProb.geometry.interpolation.nelem,1:obj.diffReacProb.geometry.interpolation.nelem,...
                    obj.diffReacProb.geometry.dvolu(:,igauss));
            end
            
            obj.coordinates = obj.diffReacProb.mesh.coord;
            obj.connectivities = obj.diffReacProb.mesh.connec;
            obj.nelem = obj.diffReacProb.geometry.interpolation.nelem;
            obj.nnode = obj.diffReacProb.geometry.interpolation.nnode;
            obj.npnod = obj.diffReacProb.geometry.interpolation.npnod;
            obj.ngaus = quadrature.ngaus;
            obj.shape = obj.diffReacProb.element.interpolation_u.shape;
        end
        
        function A_nodal_2_gauss = computeA(obj)
            A_nodal_2_gauss = sparse(obj.nelem,obj.npnod);
            fn = ones(1,obj.npnod);
            
            dirichlet_data = obj.connectivities';
            fe = zeros(obj.nnode,obj.nelem);
            
            fg = zeros(obj.ngaus,obj.nelem);
            
            for igaus = 1:obj.ngaus
                for inode = 1:obj.nnode
                    fe(inode,:) = fn(dirichlet_data(inode,:));
                    fg(igaus,:) = fg(igaus,:) + obj.shape(inode,igaus)*fe(inode,:);
                    A_nodal_2_gauss = A_nodal_2_gauss + sparse(1:obj.nelem,dirichlet_data(inode,:),ones(obj.nelem,1)*obj.shape(inode,igaus),obj.nelem,obj.npnod);
                end
            end
        end
        function P_operator=computePoperator(obj,Msmooth)
            
            dirichlet_data=zeros(obj.nnode,obj.nelem);
            for inode=1:obj.nnode
                dirichlet_data(inode,:)=obj.connectivities(:,inode);
            end
            
            T_nodal_2_gauss = sparse(obj.nelem,obj.npnod);
            
            for inode=1:obj.nnode
                T_nodal_2_gauss = T_nodal_2_gauss + sparse(1:obj.nelem,dirichlet_data(inode,:),ones(obj.nelem,1),obj.nelem,obj.npnod);
            end
            
            m = T_nodal_2_gauss*sum(Msmooth,2);
            P_operator = diag(m)\T_nodal_2_gauss;
        end
    end    
    methods (Static)
        function obj = create(settings)
            switch settings.filter
                case 'P1'
                    switch settings.optimizer
                        case {'MMA','PROJECTED GRADIENT','IPOPT'}
                            obj = Filter_P1_Density(settings.filename,settings.ptype);
                        case {'SLERP','HAMILTON-JACOBI','PROJECTED SLERP'}
                            switch settings.pdim
                                case '2D'
                                    obj = Filter_P1_LevelSet_2D(settings.filename,settings.ptype);
                                case '3D'
                                    obj = Filter_P1_LevelSet_3D(settings.filename,settings.ptype);
                            end
                    end
                case 'PDE'
                    switch settings.optimizer
                        case {'MMA','PROJECTED GRADIENT','IPOPT'}
                            obj = Filter_PDE_Density(settings.filename,settings.ptype);
                        case {'SLERP','HAMILTON-JACOBI','PROJECTED SLERP'}
                            switch settings.pdim
                                case '2D'
                                    obj = Filter_PDE_LevelSet_2D(settings.filename,settings.ptype);
                                case '3D'
                                    obj = Filter_PDE_LevelSet_3D(settings.filename,settings.ptype);
                            end
                    end
            end
        end
        
        % !! ONLY USED FOR LEVEL-SET !!
        function [F,aire] = faireF2(p,t,psi)
            np = size(p,2); nt = size(t,2);
            F = zeros(np,1);
            p1 = t(1,:); p2 = t(2,:); p3 = t(3,:);
            x1 = p(1,p1); y1 = p(2,p1); x2 = p(1,p2); y2 = p(2,p2); x3 = p(1,p3); y3 = p(2,p3);
            A = 0.5*abs((x2-x1).*(y3-y1)-(x3-x1).*(y2-y1));
            
            beta = (psi<0); 
            beta = pdeintrp(p,t,beta);
            k = find(beta>0.5);
            F = F+accumarray(p1(k)',A(k)/3',[np,1],@sum,0);
            F = F+accumarray(p2(k)',A(k)/3',[np,1],@sum,0);
            F = F+accumarray(p3(k)',A(k)/3',[np,1],@sum,0);
            aire = sum(A(k));
            
            k = find(abs(beta-1/3)<0.01);
            p1 = t(1,k); p2 = t(2,k); p3 = t(3,k);
            psi1 = psi(p1)'; psi2 = psi(p2)'; psi3 = psi(p3)';
            [psis,is] = sort([psi1;psi2;psi3],1);
            is = is+3*ones(3,1)*[0:length(k)-1];
            pl = [p1;p2;p3]; ps = pl(is);
            x1 = p(1,ps(1,:)); y1 = p(2,ps(1,:)); x2 = p(1,ps(2,:)); y2 = p(2,ps(2,:)); x3 = p(1,ps(3,:)); y3 = p(2,ps(3,:));
            x12 = (psis(1,:).*x2-psis(2,:).*x1)./(psis(1,:)-psis(2,:));
            y12 = (psis(1,:).*y2-psis(2,:).*y1)./(psis(1,:)-psis(2,:));
            x13 = (psis(1,:).*x3-psis(3,:).*x1)./(psis(1,:)-psis(3,:));
            y13 = (psis(1,:).*y3-psis(3,:).*y1)./(psis(1,:)-psis(3,:));
            A = 0.5*abs(((x12-x1).*(y13-y1)-(x13-x1).*(y12-y1)));
            F = F+accumarray(ps(1,:)',((1+psis(2,:)./(psis(2,:)-psis(1,:))+psis(3,:)./(psis(3,:)-psis(1,:))).*A/3)',[np,1],@sum,0);
            F = F+accumarray(ps(2,:)',((psis(1,:)./(psis(1,:)-psis(2,:))).*A/3)',[np,1],@sum,0);
            F = F+accumarray(ps(3,:)',((psis(1,:)./(psis(1,:)-psis(3,:))).*A/3)',[np,1],@sum,0);
            aire = aire+sum(A);
            
            k = find(abs(beta-2/3)<0.01);
            p1 = t(1,k); p2 = t(2,k); p3 = t(3,k);
            psi1 = psi(p1)'; psi2 = psi(p2)'; psi3 = psi(p3)';
            [psis,is] = sort([psi1;psi2;psi3],1,'descend');
            is = is+3*ones(3,1)*[0:length(k)-1];
            pl = [p1;p2;p3]; ps = pl(is);
            x1 = p(1,ps(1,:)); y1 = p(2,ps(1,:)); x2 = p(1,ps(2,:)); y2 = p(2,ps(2,:)); x3 = p(1,ps(3,:)); y3 = p(2,ps(3,:));
            x12 = (psis(1,:).*x2-psis(2,:).*x1)./(psis(1,:)-psis(2,:));
            y12 = (psis(1,:).*y2-psis(2,:).*y1)./(psis(1,:)-psis(2,:));
            x13 = (psis(1,:).*x3-psis(3,:).*x1)./(psis(1,:)-psis(3,:));
            y13 = (psis(1,:).*y3-psis(3,:).*y1)./(psis(1,:)-psis(3,:));
            A = 0.5*abs(((x12-x1).*(y13-y1)-(x13-x1).*(y12-y1)));
            F = F-accumarray(ps(1,:)',((1+psis(2,:)./(psis(2,:)-psis(1,:))+psis(3,:)./(psis(3,:)-psis(1,:))).*A/3)',[np,1],@sum,0);
            F = F-accumarray(ps(2,:)',((psis(1,:)./(psis(1,:)-psis(2,:))).*A/3)',[np,1],@sum,0);
            F = F-accumarray(ps(3,:)',((psis(1,:)./(psis(1,:)-psis(3,:))).*A/3)',[np,1],@sum,0);
            aire = aire-sum(A);
        end
    end
end