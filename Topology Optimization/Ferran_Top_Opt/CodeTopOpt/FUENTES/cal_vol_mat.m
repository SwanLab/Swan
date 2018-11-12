function [vol_mat,x_gp_cort,perimeter,index_border,volume,geometric_volume] = cal_vol_mat(phifunct,dim,element,problembsc,coordinates)


nelem=dim.nelem; nndof=dim.nndof; nnode=dim.nnode;
ndime=dim.ndime; npnod=dim.npnod; nunkn = dim.nunkn; nstre = dim.nstre;

node(1,:) = element.conectivities(:,1);
node(2,:) = element.conectivities(:,2);
node(3,:) = element.conectivities(:,3);
% 
sign_n(1,:) = sign(phifunct(node(1,:)));
sign_n(2,:) = sign(phifunct(node(2,:)));
sign_n(3,:) = sign(phifunct(node(3,:)));

sign_n (sign_n == 0) = 1; % we consider phi = 0 as empty

type = sum(sign_n,1);
index_empty = (type == 3);
index_border = (type == 1 | type == -1);
index_full = (type == -3);

type_border = type(index_border);
nodes_border = node(:,index_border);
phi_border  = phifunct(nodes_border);
 
pos_porder(1,:,:) = coordinates(nodes_border(1,:),:)';
pos_porder(2,:,:) = coordinates(nodes_border(2,:),:)';
pos_porder(3,:,:) = coordinates(nodes_border(3,:),:)';

[vol_void_border,x_gp_cort,perimeter] = comp_vol_case(phi_border,pos_porder,type_border);

vol_mat(1,index_empty) = 0;
vol_mat(1,index_border) = 1 - vol_void_border;
vol_mat(1,index_full) = 1;

[coordinatesn,coordinatesa] = init_coord(coordinates);
etype = element.type;
ptype = problembsc.problemtype;
[posgp,weigp,ngaus] = cal_posgp_weigp(etype,ndime,nnode,element.ngaus);
evol=zeros(ngaus,nelem);
egeometric_volum=zeros(ngaus,nelem);
%vol=zeros(npnod,1);

for igaus=1:ngaus
    [~,djacb] = cal_cartd(igaus,posgp,element,ndime,nnode,nelem,coordinatesn,coordinatesa,ptype);
    evol(igaus,:) = weigp(igaus)*(djacb.*(vol_mat'))';
    egeometric_volum(igaus,:) = weigp(igaus)*djacb';
end
geometric_volume = sum(sum(egeometric_volum));
volume = sum(sum(evol));


end



function [vol_void,x_gp_cort,perimeter] = comp_vol_case(phifunct,pos,type)

cases = [1 2; 2 3; 3 1];

for icases = 1:size(cases,1)
    index_case = sign(phifunct(cases(icases,1),:)).*sign(phifunct(cases(icases,2),:)) > 0;
    [vol_void(index_case),~,x_gp_cort(:,index_case),perimeter(index_case)] = compute_vol_void(phifunct(:,index_case),pos(:,:,index_case),type(index_case),cases(icases,1),cases(icases,2),setdiff([1:3],cases(icases,:)));
    
end

end



function [vol_void,vol_t,x_gp_cort,perimeter] = compute_vol_void(phifunct,pos,type,index_p1,index_p2,index_n)

         xcort(1,1,:) = (-phifunct(index_n,:)'.*squeeze(pos(index_p1,1,:)) + phifunct(index_p1,:)'.*squeeze(pos(index_n,1,:)))'./(phifunct(index_p1,:) - phifunct(index_n,:));
         xcort(1,2,:) = (-phifunct(index_n,:)'.*squeeze(pos(index_p1,2,:)) + phifunct(index_p1,:)'.*squeeze(pos(index_n,2,:)))'./(phifunct(index_p1,:) - phifunct(index_n,:));
         
         xcort(2,1,:) = (-phifunct(index_n,:)'.*squeeze(pos(index_p2,1,:)) + phifunct(index_p2,:)'.*squeeze(pos(index_n,1,:)))'./(phifunct(index_p2,:) - phifunct(index_n,:));
         xcort(2,2,:) = (-phifunct(index_n,:)'.*squeeze(pos(index_p2,2,:)) + phifunct(index_p2,:)'.*squeeze(pos(index_n,2,:)))'./(phifunct(index_p2,:) - phifunct(index_n,:));

         vect_1 = zeros(3,size(phifunct,2));
         vect_2 = zeros(3,size(phifunct,2));
         vect_1(1:2,:) = squeeze(xcort(1,:,:) - pos(index_n,:,:));
         vect_2(1:2,:) = squeeze(xcort(2,:,:) - pos(index_n,:,:));
         vol_h = 0.5*abs(vect_1(1,:).*vect_2(2,:)-vect_1(2,:).*vect_2(1,:));

         vectt_1 = zeros(3,size(phifunct,2));
         vectt_2 = zeros(3,size(phifunct,2));
         vectt_1(1:2,:) = pos(index_p1,:,:) - pos(index_n,:,:);
         vectt_2(1:2,:) = pos(index_p2,:,:) - pos(index_n,:,:);
         vol_t = 0.5*abs(vectt_1(1,:).*vectt_2(2,:)-vectt_1(2,:).*vectt_2(1,:));

         
         vol_void = zeros(size(vol_h));
         vol_void(1,type == 1) = (vol_t(type == 1) - vol_h(type == 1))./vol_t(type == 1); 
         vol_void(1,type == -1) = vol_h(type == -1)./vol_t(type == -1); 
         
         x_gp_cort = squeeze(mean(xcort));
         perimeter = sqrt((xcort(1,1,:) - xcort(2,1,:)).^2 + (xcort(1,2,:) - xcort(2,2,:)).^2);
end