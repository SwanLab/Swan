function structural_values = compute_structural_values(igaus,idata,structural_values,stres,strain,Ce,vstrain,dvolu,nelem,gamma,element,problembsc,dim,nstre,ftype,d_u,fextpoint)


    switch problembsc.TYPE
        
        case 'MICRO'
            tstres = structural_values.tstres;
            tstrain = structural_values.tstrain;
            matCh = structural_values.matCh;

            vol_dom = sum(dvolu(:));
            if idata ~= nstre + 1
                            
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
                switch problembsc.phisical_type
                case {'ELASTIC'}
                      matCh(istre,idata) =matCh(istre,idata) +  1/vol_dom*(stres(istre,:))*dvolu + 1/vol_dom*squeeze(Ce(istre,idata,:))'*dvolu;
                case {'THERMAL'}
                      matCh(istre,idata) =matCh(istre,idata) +  1/vol_dom*(-stres(istre,:))*dvolu + 1/vol_dom*squeeze(Ce(istre,idata,:))'*dvolu;
                end
                
            end
            
            structural_values.tstres = tstres;
            structural_values.tstrain = tstrain;
            structural_values.matCh = matCh;
            
            else
                structural_values.compliance = 0;
                structural_values.compliance_rar = 0;
                stres_m(igaus,:,:) =  stres;
                strain_m(igaus,:,:) = strain;
                for istre=1:nstre
                    for jstre=1:nstre
                        stres_m(igaus,istre,:) = squeeze(stres_m(igaus,istre,:)) + 1/vol_dom*squeeze(squeeze(Ce(istre,jstre,:))*vstrain(idata,jstre));
                    end
                    strain_m(igaus,istre,:) = shiftdim(strain_m(igaus,istre,:),2) + vstrain(idata,istre)*ones(nelem,1);
                end
                
                 for istre=1:nstre
                     
                   structural_values.compliance(igaus) = structural_values.compliance  + 1/vol_dom*(squeeze(stres_m(igaus,istre,:).*strain_m(igaus,istre,:)))'*dvolu;
                   
                   Ceinv = multinverse3x3(Ce);
                   for jstre = 1:nstre
                   structural_values.compliance_rar = structural_values.compliance_rar  + 1/vol_dom*(squeeze(stres_m(igaus,istre,:).*Ceinv(istre,jstre,:).*stres_m(igaus,jstre,:)))'*dvolu;
                   end
                   
                   structural_values.strain_hom(igaus,istre) = 1/vol_dom*(squeeze(strain_m(igaus,istre,:)))'*dvolu;
                   structural_values.stres_hom(igaus,istre) = 1/vol_dom*(squeeze(stres_m(igaus,istre,:)))'*dvolu;

                end
                
            end
  
            
            
         case 'MACRO'
            %internal work
            work1 = zeros(1,nelem);
            for istre=1:nstre
                work1 = work1 + stres(istre,:).*strain(istre,:);
            end

            switch problembsc.phisical_type
                case {'ELASTIC'}
                    switch element.material.subtype
                        case 'PLANESTRAIN'
                            trstre = stres(1,:)+stres(2,:)+stres(4,:);
                            trstra = strain(1,:)+strain(2,:);
                        case 'PLANESTRES'
                            trstre = stres(1,:)+stres(2,:);
                            trstra = strain(1,:)+strain(2,:);%+strain(4,:);
                    end
                    
                case {'THERMAL'}
                    trstre = stres(1,:)+stres(2,:);
                    trstra = strain(1,:)+strain(2,:);
            end
            
            
            
            work2 = trstre.*trstra;
            switch ftype
                case {'ELASTIC'}
                    [DtJtil] = cal_nodgfunc_objective_macro_elast(igaus,work1,work2,gamma,element,problembsc,dim,d_u,fextpoint,strain,stres);
                    structural_values.fobj(igaus,:) = work1'.*dvolu;
                    structural_values.costfunc = structural_values.costfunc + work1*dvolu;
                    
                    
                case {'THERMAL'}
                    [DtJtil] = cal_nodgfunc_objective_macro_thermal(igaus,stres,strain,gamma,element,problembsc,dim,vol_void,d_u,fextpoint);
                    structural_values.fobj(igaus,:) =  -0.5*work1'.*dvolu;
                    structural_values.costfunc = structural_values.costfunc - 0.5*work1*dvolu;
            end
            
            structural_values.DtJtil(igaus,:) = DtJtil(igaus,:);
            
            
            structural_values.stres(igaus,:,:) = stres;
            structural_values.strain(igaus,:,:) = strain;
            

    end

end


