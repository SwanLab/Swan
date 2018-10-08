classdef Filter_LevelSet_2D < Filter_LevelSet
    methods
        function obj = Filter_LevelSet_2D
            obj.max_subcells = 6;
            obj.nnodes_subelem = 3;
            obj.ndim = 2;
        end
        
        function getQuadrature_Unfitted(obj)
            obj.quadrature_unfitted = Quadrature_Triangle;
        end
        
        function getInterpolation_Unfitted(obj)
            obj.interpolation_unfitted = Triangle_Linear(obj.unfitted_mesh);
        end
        
        function createUnfittedMesh(obj)
            obj.unfitted_mesh = Mesh_Unfitted_2D(obj.mesh.duplicate,obj.diffReacProb.geometry.interpolation);
        end
        
        function M2=computeRHS_facet(obj,x,F)
%             obj.setupUnfittedMesh(x);
%             obj.unfitted_mesh.computeDvoluCut;
            
            [interp_facet,quadrature_facet] = obj.createFacet;
            interp_element = Interpolation.create(obj.mesh,obj.quadrature.order);
            
            shape_all = zeros(obj.nelem,obj.nnode);
            [~,cut_elem]=obj.findCutElements(x,obj.mesh.connec);
            
            [P_iso,active_nodes_iso]=obj.findCutPoints_Iso(x,cut_elem,obj.geometry.interpolation);
            [P_global,active_nodes_global]=obj.findCutPoints_Global(x,cut_elem);            
            
            % !! VECTORITZAR: LOOPS PETITS, ELEMENTS DIRECTES !!
%             figure, hold on
            for icut = 1:length(cut_elem)
                ielem = cut_elem(icut); inode_global = obj.mesh.connec(ielem,:);
                cutPoints_iso = P_iso(active_nodes_iso(:,:,icut),:,icut);
                cutPoints_global = P_global(active_nodes_global(:,:,icut),:,icut);
                
                connec_facets = obj.findFacetsConnectivities(cutPoints_iso,interp_element,x,inode_global);
                
                for i = 1:size(connec_facets,1)
                    for igaus = 1:quadrature_facet.ngaus
                        for idime = 1:interp_element.ndime
                            facet_posgp(igaus,idime) = interp_facet.shape(igaus,:)*cutPoints_iso(connec_facets(i,:),idime);
                        end
                    end
                    interp_element.computeShapeDeriv(facet_posgp');
                    facet_deriv(:,:) = interp_facet.deriv(:,:,:);
                    
                    % !! How mapping is done for 2D cases??? !!
                    t = [0; norm(diff(cutPoints_global(connec_facets(i,:),:)))];
                    dt_dxi = (facet_deriv'*t)/interp_facet.dvolu;
                    
                    f = (interp_element.shape*quadrature_facet.weigp')'*F(inode_global)/interp_facet.dvolu;
                    shape_all(ielem,:) = shape_all(ielem,:) + (interp_element.shape*(dt_dxi.*quadrature_facet.weigp')*f)';
                    
%                     plot(obj.coordinates(obj.mesh.connec(ielem,:),1),obj.mesh.coord(obj.mesh.connec(ielem,:),2),'.-b'); plot(obj.mesh.coord(obj.mesh.connec(ielem,[1 3]),1),obj.mesh.coord(obj.mesh.connec(ielem,[1 3]),2),'.-b');
%                     plot(cutPoints_global(connec_facets(i,:),1),cutPoints_global(connec_facets(i,:),2),'-xr');
%                     title('Cut Elements & Cut points in GLOBAL coordinates'), axis('equal')
                end
            end
            M2=obj.rearrangeOutputRHS(shape_all);
        end
        
        function [P,active_nodes]=findCutPoints_Iso(obj,x,cut_elem,interpolation)
            gamma_1=permute(x(obj.mesh.connec(cut_elem,:)),[2 3 1]);
            gamma_2=permute([x(obj.mesh.connec(cut_elem,2:end)),x(obj.mesh.connec(cut_elem,1))],[2 3 1]);
            P1=repmat(interpolation.pos_nodes,[1 1 size(cut_elem)]);
            P2=repmat([interpolation.pos_nodes(2:end,:);interpolation.pos_nodes(1,:)],[1 1 size(cut_elem)]);
            P=P1+gamma_1.*(P2-P1)./(gamma_1-gamma_2);
            active_nodes = sign(gamma_1.*gamma_2)<0;
        end
        
        function [P,active_nodes]=findCutPoints_Global(obj,x,cut_elem,interpolation)
            index1 = permute(obj.mesh.connec(cut_elem,:),[2 3 1]);
            index2 = [permute(obj.mesh.connec(cut_elem,2:end),[2 3 1]);...
                permute(obj.mesh.connec(cut_elem,1),[2 3 1])];
            gamma_1=x(index1);
            gamma_2=x(index2);
            coord1 = obj.mesh.coord(:,1); coord2 = obj.mesh.coord(:,2);
            P1=[coord1(index1) coord2(index1)];
            P2=[coord1(index2) coord2(index2)];
            P=P1+gamma_1.*(P2-P1)./(gamma_1-gamma_2);
            active_nodes = sign(gamma_1.*gamma_2)<0;
        end
        
        function [interp_facet,quadrature_facet] = createFacet(obj)
            quadrature_facet = Quadrature.set('LINE');
            interp_facet = Line_Linear;
            quadrature_facet.computeQuadrature(obj.quadrature.order);
            interp_facet.computeShapeDeriv(quadrature_facet.posgp);
        end
    end
    
    methods (Static)
        function A=computeDvoluCut(elcrd)
            x1 = elcrd(:,1,1); y1 = elcrd(:,1,2); x2 = elcrd(:,2,1); y2 = elcrd(:,2,2); x3 = elcrd(:,3,1); y3 = elcrd(:,3,2);
            A = 0.5*abs((x2-x1).*(y3-y1)-(x3-x1).*(y2-y1));
        end
        
        function facets_connectivities = findFacetsConnectivities(elem_cutPoints_iso,interpolation,x,inode_global)
            if size(elem_cutPoints_iso,1) == 2
                facets_connectivities = [1 2];
            elseif size(elem_cutPoints_iso,1) == 4
                del_coord = [interpolation.pos_nodes;elem_cutPoints_iso];
                DT=delaunayTriangulation(del_coord);
                del_connec=DT.ConnectivityList;
                
                node_positive_iso = find(x(inode_global)>0);
                
                for idel = 1:length(node_positive_iso)
                    [a, ~] = find(del_connec==node_positive_iso(idel));
                    facets_connectivities(idel,:) = del_connec(a(end),del_connec(a(end),:)~=node_positive_iso(idel))-interpolation.nnode;
                end
            else
                error('Case still not implemented.')
            end
        end
        
        function djacob = mapping(elem_cutPoints_global,facets_connectivities,facet_deriv,dvolu)
            t = [0; norm(diff(elem_cutPoints_global(facets_connectivities,:)))];
            djacob = (facet_deriv'*t)/dvolu;
        end
    end
end

