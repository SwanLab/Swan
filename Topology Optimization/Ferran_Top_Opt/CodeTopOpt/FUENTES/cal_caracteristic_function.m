function [chi] = cal_caracteristic_function(igaus,phifunct,element,dim,problembsc,vol_void)
% Compute caracteristic function chi, for each element. 
% chi is defined in terms of the levet set function.
% If phi(x)<0 chi=1; if phi(x)>0 chi=epsi (where epsi <<1)
%

nelem=dim.nelem;  nnode=dim.nnode; ndime=dim.ndime; npnod= dim.npnod;
epsi = element.material.opt_epsi;
if (problembsc.smoothChi==0)
    % Chi se define en cada pto de gauss via phigp(x_kg) 
    [phigp] = interpol(phifunct,element,dim,problembsc);
    chi = ones(1,nelem);
    outside= (phigp(igaus,:) > 0);
    if any(outside)
        chi(outside)= epsi;
    end

    chi(element.boundary_elements) = 1;
    
elseif (problembsc.smoothChi==1)
    % Chi se interpona a partir de un chi nodal que se define via phi_nodal
    nod_chi = ones(npnod,1);
    outside = (phifunct > 0);
    if any(outside)
        nod_chi(outside)= epsi;
    end
    
    nod_chi(unique(element.conectivities(element.boundary_elements,:))) = 1;
        
    [chigp] = interpol(nod_chi,element,dim,problembsc);
    chi = chigp(igaus,:);
     
elseif (problembsc.smoothChi==2)
    
    chi = 1-vol_void;
   % chi(abs(chi)<1e-12) =  element.material.opt_epsi;
    %chi(:,element.boundary_elements) = 1;
    
elseif (problembsc.smoothChi==3)
    dirichlet_data=zeros(nnode,nelem);
    ephi=zeros(nnode,nelem);
    chi = ones(1,nelem);
    for i=1:nnode
        dirichlet_data(i,:)= element.conectivities(:,i);
        ephi(i,:)= phifunct(dirichlet_data(i,:));
    end
    avgphi = mean(ephi);
    outside= (avgphi > 0);
    if any(outside)
        chi(outside)= epsi;
    end
    chi(element.boundary_elements) = 1;
end





end

