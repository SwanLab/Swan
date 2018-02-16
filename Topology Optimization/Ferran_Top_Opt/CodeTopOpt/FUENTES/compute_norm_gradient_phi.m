function norm_grad_phi = compute_norm_gradient_phi(phi,element,dim,problembsc,coordinates,index_border)

        nelem=dim.nelem;  nnode=dim.nnode;
        ndime = dim.ndime; nstre = dim.ndime; nunkn = 1;
        ptype = problembsc.problemtype;
        ftype = 'THERMAL';
        etype = element.type;
        [posgp] = cal_posgp_weigp(element.type,ndime,nnode,element.ngaus);

        for igaus=1:element.ngaus
        [cartd,~] = cal_cartd(igaus,posgp,element,ndime,nnode,nelem,coordinates,coordinates,ptype);
        [Bmat] = cal_B(cartd,nstre,nnode,nunkn,nelem,etype,ptype,ftype);
        end

        grad_phi = zeros(sum(index_border),ndime);
        for inode = 1:dim.nnode
            lnodes = element.conectivities(index_border,inode);
            for idime = 1:dim.ndime*nunkn
                grad_phi(:,idime) = grad_phi(:,idime) + squeeze(Bmat(idime,inode,index_border)).*phi(lnodes);
            end
        end
        norm_grad_phi(:,1) = sqrt(grad_phi(:,1).^2 + grad_phi(:,2).^2);




end