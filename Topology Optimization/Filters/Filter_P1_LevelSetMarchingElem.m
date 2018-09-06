classdef Filter_P1_LevelSetMarching < Filter_P1
    properties
        quadrature
        geometry        
        quadrature_del
        interp_del
        shape_full
    end
    methods    
        function obj = Filter_P1_LevelSetMarching(problemID,scale)
            obj@Filter_P1(problemID,scale);
        end
        function preProcess(obj)
            preProcess@Filter_P1(obj)
            obj.quadrature = Quadrature.set(obj.diffReacProb.geometry.type);
            mesh=obj.diffReacProb.mesh.duplicate;
            obj.geometry= Geometry(mesh,'LINEAR');    
            mesh_del=mesh.duplicate; 
            switch mesh.pdim
                case '2D'
                    mesh_del.geometryType='TRIANGLE';
                    obj.quadrature_del=Quadrature_Triangle;
                    obj.quadrature_del.computeQuadrature('QUADRATIC');
                    obj.interp_del=Triangle_Linear(mesh_del);
                    obj.interp_del.computeShapeDeriv(obj.quadrature_del.posgp)
                case '3D'
                    mesh_del.geometryType='TETRAHEDRA';
                    obj.quadrature_del=Quadrature_Tetrahedra;
                    obj.quadrature_del.computeQuadrature('LINEAR');
                    obj.interp_del=Tetrahedra(mesh_del);
                    obj.interp_del.computeShapeDeriv(obj.quadrature_del.posgp)                    
            end   
            obj.initGeometry
            obj.shape_full=obj.integrateFull;
        end
        function x_gp = getP0fromP1(obj,x)
            if norm(x) == norm(obj.x)
                x_gp=obj.x_reg;
            else
                switch obj.geometry.type
                    case 'TRIANGLEd'
                        M2=obj.faireF2(obj.coordinates',obj.connectivities',x);
                    otherwise
                        % tic
                        % M2_faire=obj.faireF2(obj.coordinates',obj.connectivities',x);
                        %  toc
                        % tic
                        M2=obj.computeRHS(x);
                        %   toc
                        %error=abs(norm(M2)-norm(M2_faire))/norm(M2_faire)
                end
                
                x_gp = obj.P_operator*M2;
                obj.x_reg=x_gp;
            end
        end
        function initGeometry(obj)
            obj.quadrature.computeQuadrature('LINEAR');
            obj.geometry.interpolation.computeShapeDeriv(obj.quadrature.posgp)
            obj.geometry.computeGeometry(obj.quadrature,obj.geometry.interpolation);
        end
        function [full_elem,cut_elem]=findCutElements(obj,x)
            phi_nodes=x(obj.connectivities);
            phi_case=sum((sign(phi_nodes)<0),2);
            
            full_elem = phi_case==size(obj.connectivities,2);
            null_elem = phi_case==0;
            indexes = (1:size(obj.connectivities,1))';
            cut_elem = indexes(~(full_elem+null_elem));
        end
        function shape=integrateFull(obj)
            shape=zeros(size(obj.connectivities,1),size(obj.connectivities,2));
            for igauss=1:size(obj.geometry.interpolation.shape,2)
                shape=shape+obj.geometry.interpolation.shape(:,igauss)'.*obj.geometry.dvolu(:,igauss);
            end
        end
        function [P,active_nodes]=findCutPoints(obj,x,cut_elem)
            switch obj.diffReacProb.mesh.pdim
                case '2D'
                    gamma_1=permute(x(obj.connectivities(cut_elem,:)),[2 3 1]);
                    gamma_2=permute([x(obj.connectivities(cut_elem,2:end)),x(obj.connectivities(cut_elem,1))],[2 3 1]);
                    P1=repmat(obj.geometry.interpolation.pos_nodes,[1 1 size(cut_elem)]);
                    P2=repmat([obj.geometry.interpolation.pos_nodes(2:end,:);obj.geometry.interpolation.pos_nodes(1,:)],[1 1 size(cut_elem)]);
                    P=P1+gamma_1.*(P2-P1)./(gamma_1-gamma_2);
                    active_nodes = sign(gamma_1.*gamma_2)<0;
%                     gamma_1=x(obj.connectivities(cut_elem,:));
%                     gamma_2=[x(obj.connectivities(cut_elem,2:end)),x(obj.connectivities(cut_elem,1))];
%                     active_nodes = sign(gamma_1.*gamma_2)<0;
%                     P1=repmat(obj.geometry.interpolation.pos_nodes,[1 1 size(cut_elem)]);
%                     P2=repmat([obj.geometry.interpolation.pos_nodes(2:end,:);obj.geometry.interpolation.pos_nodes(1,:)],[1 1 size(cut_elem)]);
%                     P=P1+gamma_1.*(P2-P1)./(gamma_1-gamma_2);
                    
                case '3D'
                   % iteration_1=[1 1 1 2 2 3];
                   % iteration_2=[2 3 4 3 4 4];
                   iteration_1=obj.geometry.interpolation.iteration(1,:);
                   iteration_2=obj.geometry.interpolation.iteration(2,:);
                    gamma_1=permute(x(obj.connectivities(cut_elem,iteration_1)),[2 3 1]);
                    gamma_2=permute(x(obj.connectivities(cut_elem,iteration_2)),[2 3 1]);
                    P1=repmat(obj.geometry.interpolation.pos_nodes(iteration_1,:),[1 1 size(cut_elem)]);
                    P2=repmat(obj.geometry.interpolation.pos_nodes(iteration_2,:),[1 1 size(cut_elem)]);
                    P=P1+gamma_1.*(P2-P1)./(gamma_1-gamma_2);
                    active_nodes = sign(gamma_1.*gamma_2)<0;
            end
        end
        function A=computeDvoluCut(obj,elcrd)
            switch obj.diffReacProb.mesh.pdim
                case '2D'
                    x1 = elcrd(:,1,1); y1 = elcrd(:,1,2); x2 = elcrd(:,2,1); y2 = elcrd(:,2,2); x3 = elcrd(:,3,1); y3 = elcrd(:,3,2);
                    A = 0.5*abs((x2-x1).*(y3-y1)-(x3-x1).*(y2-y1));
                case '3D'
                    x1 = elcrd(:,1,1); y1 = elcrd(:,1,2); z1=elcrd(:,1,3);
                    x2 = elcrd(:,2,1); y2 = elcrd(:,2,2); z2=elcrd(:,2,3);
                    x3 = elcrd(:,3,1); y3 = elcrd(:,3,2); z3=elcrd(:,3,3);
                    x4 = elcrd(:,4,1); y4 = elcrd(:,4,2); z4=elcrd(:,4,3);
                    J=x1.*y3.*z2-x1.*y2.*z3+x2.*y1.*z3-x2.*y3.*z1-x3.*y1.*z2+x3.*y2.*z1+x1.*y2.*z4-x1.*y4.*z2-x2.*y1.*z4+x2.*y4.*z1+...
                    x4.*y1.*z2-x4.*y2.*z1-x1.*y3.*z4+x1.*y4.*z3+x3.*y1.*z4-x3.*y4.*z1-x4.*y1.*z3+x4.*y3.*z1+x2.*y3.*z4-x2.*y4.*z3...
                    -x3.*y2.*z4+x3.*y4.*z2+x4.*y2.*z3-x4.*y3.*z2;
                    A=J/6;
            end
        end
        function shape_all=integrateCut(obj,phi_cut, global_connec, A, shape_all,pos_gp_del_natural)
            notcompute = phi_cut > 0;
            isnegative = ~any(notcompute');
            v=zeros(size(global_connec,1),size(shape_all,2));
             
            dvolu=sum(obj.geometry.dvolu,2)/obj.geometry.interpolation.dvolu;
            for igauss=1:obj.quadrature_del.ngaus
                obj.geometry.interpolation.computeShapeDeriv(pos_gp_del_natural(:,:,igauss)');
                v=v+isnegative'.*(obj.geometry.interpolation.shape'.*A.*dvolu(global_connec))*obj.quadrature_del.weigp(igauss);
            end
            for idelaunay=1:size(v,2)
                shape_all(:,idelaunay)=shape_all(:,idelaunay)+accumarray(global_connec,v(:,idelaunay),[obj.nelem,1],@sum,0);
            end
        end
        function pos_gp_del_natural=computePosGpDelaunayNatural(obj,elcrd)
            pos_gp_del_natural=zeros(size(elcrd,1),size(elcrd,3));
            for igauss=1:size(obj.interp_del.shape,2)
                for idime=1:size(elcrd,3)
                    pos_gp_del_natural(:,idime,igauss)=elcrd(:,:,idime)*obj.interp_del.shape(:,igauss);
                end
            end
        end
        function shape_all=computeFullElements(obj,full_elem)
            shape_all=zeros(size(obj.connectivities,1),size(obj.connectivities,2));
            shape_all(full_elem,:)=obj.shape_full(full_elem,:);
        end
        function [elecoord,global_connec,phi_cut]=MarchingCubes(obj,x,cut_elem)
                [P,active_nodes]=obj.findCutPoints(x,cut_elem);
                elecoord=[];phi_cut=[];global_connec=[];xcoord=[];ycoord=[];
                coord_cut=repmat(obj.geometry.interpolation.pos_nodes,[1 1 size(cut_elem)]);
                P=permute(P,[1 3 2]);
                active_nodes=permute(active_nodes,[1 3 2]);
                coord_cut=permute(coord_cut,[1 3 2]);
                a=obj.geometry.interpolation.nnode;
                for inode=1:size(P,1)
                    for idime=1:obj.geometry.interpolation.ndime
                        coord_cut(a+inode,active_nodes(inode,:),idime)=P(inode,active_nodes(inode,:),idime);
                    end
                end
                    
        %             cutnodes=P(active_nodes(:,idime,:)
                   %  coord_cut(4,:,:) = [repmat(obj.geometry.interpolation.pos_nodes,[1 1 size(cut_elem)]);P(active_nodes(:,:,ielem),:,ielem)];
                   
            %    end

            num=repmat([1:obj.geometry.interpolation.nnode],[size(cut_elem) 1]);
            cutcases=x(obj.connectivities(cut_elem,:))<0;
            test1=num.*cutcases;
            test2=sum(test1,2);
            test3=sum(cutcases,2);
            cases_list=diag(obj.geometry.interpolation.selectcases(test2,test3));
            connec=obj.geometry.interpolation.cases(cases_list);
            
            for ielem=1:size(cut_elem)
                triangles=connec{ielem};
                x_elem=[x(obj.connectivities(cut_elem(ielem),:));zeros(sum(active_nodes(:,ielem)),1)];
                triangles_coordx=  [obj.geometry.interpolation.pos_nodes(:,1);P(active_nodes(:,ielem),ielem,1)];
                triangles_coordy=  [obj.geometry.interpolation.pos_nodes(:,2);P(active_nodes(:,ielem),ielem,2)];
               % triangles_coordx=coord_cut(:,ielem,1);
               % triangles_coordy=coord_cut(:,ielem,2);
%                if ielem<5
              %  figure
             %    triplot(triangles,triangles_coordx,triangles_coordy)
%                end
                xcoord=[xcoord;triangles_coordx(triangles)];
                ycoord=[ycoord;triangles_coordy(triangles)];
                phi_cut =[phi_cut;x_elem(triangles)];
                global_connec=[global_connec;repmat(cut_elem(ielem),[size(triangles,1) 1])];
            end
            elecoord(:,:,1)=xcoord;
            elecoord(:,:,2)=ycoord;
%                     del_x=[x(obj.connectivities(cut_elem(ielem),:));zeros(size(P(active_nodes(:,:,ielem)),1),1)];
%                     DT=delaunayTriangulation(del_coord);
%                     del_connec=DT.ConnectivityList;
%                     for idelaunay=1:size(del_connec,1)
%                         elecoord=[elecoord;permute(del_coord(del_connec(idelaunay,:),:)',[3 2 1])];
%                         
%                     end
%                     
                
        end
        function [elecoord,global_connec,phi_cut]=computeDelaunay(obj,x,cut_elem)
            [P,active_nodes]=obj.findCutPoints(x,cut_elem);
            elecoord=[];phi_cut=[];global_connec=[];
            for ielem=1:length(cut_elem)
                del_coord = [obj.geometry.interpolation.pos_nodes;P(active_nodes(:,:,ielem),:,ielem)];
                del_x=[x(obj.connectivities(cut_elem(ielem),:));zeros(size(P(active_nodes(:,:,ielem)),1),1)];
                DT=delaunayTriangulation(del_coord);
                del_connec=DT.ConnectivityList;

             %   figure(4)
             %   triplot(DT)

                for idelaunay=1:size(del_connec,1)
                    elecoord=[elecoord;permute(del_coord(del_connec(idelaunay,:),:)',[3 2 1])];
                    phi_cut =[phi_cut;del_x(del_connec(idelaunay,:))'];
                end
                global_connec=[global_connec;repmat(cut_elem(ielem),[size(del_connec,1) 1])];
            end
        end
        function M2=rearrangeOutputRHS(obj,shape_all)
              M2=zeros(obj.npnod,1);
            for inode=1:obj.nnode
                p = obj.connectivities(:,inode);
                M2 = M2+accumarray(p,shape_all(:,inode),[obj.npnod,1],@sum,0);
            end
        end
        function M2=computeRHS(obj,x)
            [full_elem,cut_elem]=obj.findCutElements(x);
            shape_all=obj.computeFullElements(full_elem);
            shape_all2=shape_all;
            if ~isempty(cut_elem)                
                [delaunaycoord,global_connec,phi_cut]=obj.MarchingCubes(x,cut_elem);
                [delaunaycoord2,global_connec2,phi_cut2]=obj.computeDelaunay(x,cut_elem);
                dvolu_cut=obj.computeDvoluCut(delaunaycoord);
                pos_gp_del_natural=obj.computePosGpDelaunayNatural(delaunaycoord);                                             
                shape_all=obj.integrateCut(phi_cut, global_connec, dvolu_cut, shape_all,pos_gp_del_natural);
                
                dvolu_cut2=obj.computeDvoluCut(delaunaycoord2);
                pos_gp_del_natural2=obj.computePosGpDelaunayNatural(delaunaycoord2);                             
                shape_all2=obj.integrateCut(phi_cut2, global_connec2, dvolu_cut2, shape_all2,pos_gp_del_natural2);
            end
            M2=obj.rearrangeOutputRHS(shape_all);   
            M22=obj.rearrangeOutputRHS(shape_all2);
            error=abs(norm(M2)-norm(M22))/norm(M2)
        end
    end
end

