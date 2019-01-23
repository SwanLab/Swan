function [DtC] = cal_DtC(dim,gamma_star,tstres,tstrain,ngaus,phigp,element,Msmooth,coordinates,problembsc,phifunct,vol_void,d_u,fextpoint)
% gamma = element.material.opt_epsi


[coordinatesn,coordinatesa] = init_coord(coordinates);


xg = interpol(coordinates(:,1),element,dim,problembsc);
yg = interpol(coordinates(:,2),element,dim,problembsc);

nelem=dim.nelem; nndof=dim.nndof; nnode=dim.nnode;
ndime=dim.ndime; npnod=dim.npnod; nunkn = dim.nunkn; nstre = dim.nstre;
ftype = problembsc.phisical_type;
ptype = problembsc.problemtype;

[posgp,weigp,ngaus] = cal_posgp_weigp(element.type,ndime,nnode,element.ngaus);


[h1_menos,h2_menos] = coef_DtC(gamma_star,element,problembsc);
[h1_mas,h2_mas] = coef_DtC(1/gamma_star,element,problembsc);
%h1_menos = -h1_menos; h2_menos = -h2_menos;
DtC = zeros(nstre,nstre,ngaus,nelem);

% I = [1 3; 3 2];
% n=2;
% for i=1:n
%     for j=1:n
%         for k=1:n
%             for l=1:n
%                 a=I(i,j);
%                 b=I(k,l);
%                 astre=zeros(ngaus,nelem);
%                 bstre=zeros(ngaus,nelem);
%                 nsigma=zeros(ngaus,nelem);
%                 atrace=zeros(ngaus,nelem);
%                 btrace=zeros(ngaus,nelem);
%                 
%                 for istre=1:nstre
%                     astre(:,:)=tstres(a,:,istre,:);
%                     bstre(:,:)=tstres(b,:,istre,:);
%                     if (istre==3)
%                         nsigma = nsigma + 2*astre.*bstre;
%                     else
%                         nsigma = nsigma + astre.*bstre;
%                     end
%                     if (istre~=3)
%                         atrace = atrace + astre;
%                         btrace = btrace + bstre;
%                     end
%                 end
%                 
%                 DtC_menos = h1_menos.*nsigma + h2_menos.*atrace.*btrace;
%                 DtC_mas = h1_mas.*nsigma + h2_mas.*atrace.*btrace;
%                 DtC(a,b,:,:) = vol_void.*(DtC_mas) + (1-vol_void).*(-DtC_menos);
% 
%             end
%         end
%     end
% end

switch problembsc.phisical_type
    case {'ELASTIC'}
        
        for i=1:nstre
            for j=1:nstre
                nsigma=zeros(ngaus,nelem);
                itrace = zeros(ngaus,nelem);
                jtrace = zeros(ngaus,nelem);
                Voigt_contr_tensor_prod_coef = [1 1 2];
                Voigt_trace_coef = [1 1 0];
                for istre=1:nstre
                    
                    nsigma = nsigma + shiftdim(Voigt_contr_tensor_prod_coef(istre)*tstres(i,:,istre,:).*tstres(j,:,istre,:),2);
                    itrace = itrace + shiftdim(Voigt_trace_coef(istre)*tstres(i,:,istre,:),2);
                    jtrace = jtrace + shiftdim(Voigt_trace_coef(istre)*tstres(j,:,istre,:),2);
                end
                
                DtC_menos = h1_menos.*nsigma + h2_menos.*itrace.*jtrace;
                DtC_mas = h1_mas.*nsigma + h2_mas.*itrace.*jtrace;
                DtC(i,j,:,:) = vol_void.*(DtC_mas) + (1-vol_void).*(-DtC_menos);
                
            end
        end
        
        
    case {'THERMAL'}
        
        
        for i=1:nstre
            for j=1:nstre
                nsigma=zeros(ngaus,nelem);
                Voigt_contr_tensor_prod_coef = [1 1];
                for istre=1:nstre
                    nsigma = nsigma + shiftdim(Voigt_contr_tensor_prod_coef(istre)*tstrain(i,:,istre,:).*tstrain(j,:,istre,:),2);
                   
                end
                
               % DtC_menos = h1_menos.*nsigma;
               % DtC_mas = h1_mas.*nsigma;
                
                DtC_menos = (1-vol_void).*h1_menos.*nsigma;
                DtC_mas = (1-vol_void).*h1_mas.*nsigma;
                
                DtC(i,j,:,:) = vol_void.*(DtC_mas) + (1-vol_void).*(-DtC_menos);
                
            end
        end
        
        
end


% for i=1:nstre
%     for j=1:nstre
%         
%         for igaus = 1:ngaus
%                 [cartd,djacb] = cal_cartd(igaus,posgp,element,ndime,nnode,nelem,coordinatesn,coordinatesa,ptype);
%                     dvolu = weigp(igaus)*djacb;
% 
%             stres = shiftdim(tstres(i,igaus,:,:),2);
%             strain = shiftdim(tstrain(j,igaus,:,:),2);
%             
%             work1 = zeros(1,nelem);
%             for istre=1:nstre
%                 work1 = work1 + stres(istre,:).*strain(istre,:);
%             end
%             
%             switch problembsc.phisical_type
%                 case {'ELASTIC'}
%                     switch element.material.subtype
%                         case 'PLANESTRAIN'
%                             trstre = stres(1,:)+stres(2,:)+stres(4,:);
%                             trstra = strain(1,:)+strain(2,:);
%                         case 'PLANESTRES'
%                             trstre = stres(1,:)+stres(2,:);
%                             trstra = strain(1,:)+strain(2,:)+strain(4,:);
%                     end
%                     
%                 case {'THERMAL'}
%                     trstre = stres(1,:)+stres(2,:);
%                     trstra = strain(1,:)+strain(2,:);
%             end
%             
%             work2 = trstre.*trstra;
%             switch ftype
%                 case {'ELASTIC'}
%                     [DtJtil] = cal_nodgfunc_objective_macro_elast(igaus,work1,work2,phifunct,element,problembsc,dim,vol_void,d_u,fextpoint);
%  %                   structural_values.fobj(igaus,:) = work1'.*dvolu;
%  %                   structural_values.costfunc = structural_values.costfunc + work1*dvolu;
%                 case {'THERMAL'}
%                     [DtJtil] = cal_nodgfunc_objective_macro_thermal(igaus,stres,strain,phifunct,element,problembsc,dim,vol_void,d_u,fextpoint);
%  %                   structural_values.fobj(igaus,:) = -0.5*work1'.*dvolu;
%   %                  structural_values.costfunc = structural_values.costfunc - 0.5*work1*dvolu;
%             end
%             
%             structural_values.DtJtil(i,j,igaus,:) = DtJtil(igaus,:);
%             
%         end
%         
%     end
% end



end



