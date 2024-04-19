function [volume,vol_void] = cal_omega(phifunct,dim,element,problembsc,coordinates)


nelem=dim.nelem; nndof=dim.nndof; nnode=dim.nnode;
ndime=dim.ndime; npnod=dim.npnod; nunkn = dim.nunkn; nstre = dim.nstre;

node(1,:) = element.conectivities(:,1);
node(2,:) = element.conectivities(:,2);
node(3,:) = element.conectivities(:,3);

sign_n(1,:) = sign(phifunct(node(1,:))) + (phifunct(node(1,:))==0);
sign_n(2,:) = sign(phifunct(node(2,:))) + (phifunct(node(2,:))==0);
sign_n(3,:) = sign(phifunct(node(3,:))) + (phifunct(node(3,:))==0);

type = sum(sign_n,1);
index_empty = type == 3;
index_border = (type == 1 | type == -1);
index_full = type == -3;

type_border = type(index_border);
nodes_border = node(:,index_border);
phi_border  = phifunct(nodes_border);

pos_porder(1,:,:) = coordinates(nodes_border(1,:),:)';
pos_porder(2,:,:) = coordinates(nodes_border(2,:),:)';
pos_porder(3,:,:) = coordinates(nodes_border(3,:),:)';

vol_void = zeros(1,nelem);
vol_void(1,index_empty) = 1;
vol_void(1,index_border) = comp_vol_case(phi_border,pos_porder,type_border);
vol_void(1,index_full) = 0;
vol_void(1,element.boundary_elements) = 0; %fixed elements

%[coordinatesn,coordinatesa] = init_coord(coordinates);
%etype = element.type;
%ptype = problembsc.problemtype;
%[posgp,weigp,ngaus] = cal_posgp_weigp(etype,ndime,nnode,element.ngaus);
%evol=zeros(ngaus,nelem);
%geometric_volum=zeros(ngaus,nelem);
%vol=zeros(npnod,1);

area = pdetrg(coordinates',element.conectivities');

% for igaus=1:ngaus
%     [~,djacb] = cal_cartd(igaus,posgp,element,ndime,nnode,nelem,coordinatesn,coordinatesa,ptype);
%     evol(igaus,:) = weigp(igaus)*(djacb.*(1-vol_void'))';
%     geometric_volum(igaus,:) = weigp(igaus)*djacb';
% end

%volume = sum(sum(evol))/sum(sum(geometric_volum));
evol = area.*(1-vol_void);
volume = sum(evol)/sum(area);


end



function vol_void = comp_vol_case(phifunct,pos,type)

cases = [1 2; 2 3; 3 1];

for icases = 1:size(cases,1)

    index_case = ( sign(phifunct(cases(icases,1),:)) + (phifunct(cases(icases,1),:)==0) ).* ( sign(phifunct(cases(icases,2),:)) + (phifunct(cases(icases,2),:)==0))> 0;
    vol_void(index_case) = compute_vol_void(phifunct(:,index_case),pos(:,:,index_case),type(index_case),cases(icases,1),cases(icases,2),setdiff([1:3],cases(icases,:)));
end

end



function [vol_void,vol_t] = compute_vol_void(phifunct,pos,type,index_p1,index_p2,index_n)
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
 
end
