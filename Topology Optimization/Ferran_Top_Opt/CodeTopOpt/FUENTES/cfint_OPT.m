function [force_int,structural_values] = cfint_OPT(dim,...
    coordinatesn,coordinatesa,element,problembsc,gamma,idata,structural_values,vstrain,d_u,fextpoint)

nelem=dim.nelem; nndof=dim.nndof; nnode=dim.nnode;
ndime=dim.ndime; npnod=dim.npnod; nunkn = dim.nunkn; nstre = dim.nstre;

dirichlet_data = zeros(nnode,nelem);
idx   = zeros(nnode*nunkn,nelem);
force_int  = zeros(nndof,1);
ptype = problembsc.problemtype;
ftype = problembsc.phisical_type;

for i=1:nnode
    dirichlet_data(i,:)= element.conectivities(:,i);
end

etype = element.type;
[posgp,weigp,ngaus] = cal_posgp_weigp(element.type,ndime,nnode,element.ngaus);


for inode=1:nnode
    for idime=1:nunkn
        idx(nunkn*inode-nunkn+idime,:) = nunkn.*dirichlet_data(inode,:)-nunkn+idime;
    end
end


for igaus=1:ngaus
    efint = zeros(nnode*nunkn,nelem);
    [cartd,djacb] = cal_cartd(igaus,posgp,element,ndime,nnode,nelem,coordinatesn,coordinatesa,ptype);
    [Bmat] = cal_B(cartd,nstre,nnode,nunkn,nelem,etype,ptype,ftype);
    dvolu = weigp(igaus)*djacb;
    
    % [chi] = cal_caracteristic_function(igaus,gamma,element,dim,problembsc,vol_void);
    %chi = 1-vol_void;
    %chi(abs(chi)<1e-12) =  element.material.opt_epsi;
     [Ce] = compute_consitutive_law(element,problembsc,igaus,dim,gamma);
     [stres,strain] = compute_stres_strain(dim,Bmat,element,problembsc,coordinatesn,coordinatesa,Ce,d_u);
     

     
     
     
    
    switch element.material.subtype
        case {'PLANESTRAIN','PLANESTRES'}
            for iv=1:nnode*nunkn
                ivBmat(1:nstre,1:nelem) = Bmat(1:nstre,iv,1:nelem);
                for istre=1:nstre
                    efint(iv,:)=efint(iv,:)+ivBmat(istre,:).*stres(istre,:).*dvolu';
                end
            end
        case 'AXI'
    end
    % assembling of the internal forces
    for k=1:nunkn*nnode
        %size(sparse(idx(k,:),1,efint(k,:),nndof,1));
        force_int = force_int + sparse(idx(k,:),1,efint(k,:),nndof,1);
    end
    
    switch problembsc.TYPE
        
        case 'MICRO'
            tstres = structural_values.tstres;
            tstrain = structural_values.tstrain;
            matCh = structural_values.matCh;

            vol_dom = sum(dvolu(:));
            
                            
            tstres(idata,igaus,:,:) = stres;
            tstrain(idata,igaus,:,:) = strain;

                for istre=1:nstre
                  for jstre=1:nstre
                    tstres(idata,igaus,istre,:) = squeeze(tstres(idata,igaus,istre,:)) + 1/vol_dom*squeeze(squeeze(Ce(istre,jstre,:))*vstrain(idata,jstre));
                  end  
                  tstrain(idata,igaus,istre,:) = shiftdim(tstrain(idata,igaus,istre,:),3) + vstrain(idata,istre)*ones(nelem,1);
                end


          % contribucion a la C homogeneizada
            for istre=1:nstre
                switch ftype
                    case {'ELASTIC'}
                      matCh(istre,idata) =matCh(istre,idata) +  1/vol_dom*(stres(istre,:))*dvolu + 1/vol_dom*squeeze(Ce(istre,idata,:))'*dvolu;
                    case {'THERMAL'}
                      matCh(istre,idata) =matCh(istre,idata) +  1/vol_dom*(-stres(istre,:))*dvolu + 1/vol_dom*squeeze(Ce(istre,idata,:))'*dvolu;
                end

            end

            structural_values.tstres = tstres;
            structural_values.tstrain = tstrain;
            structural_values.matCh = matCh;

            
        case 'MACRO'
            structural_values.stres(igaus,:,:) = stres;
            structural_values.strain(igaus,:,:) = strain;
            
    end
    
    structural_values.mu = element.mu_func(gamma);                                             
    structural_values.kappa = element.kappa_func(gamma);  
end



end