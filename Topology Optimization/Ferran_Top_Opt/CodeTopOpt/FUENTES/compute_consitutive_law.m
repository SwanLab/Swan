function [Ce] = compute_consitutive_law(element,problembsc,igaus,dim,gamma)    

nelem=dim.nelem; nndof=dim.nndof; nnode=dim.nnode;
ndime=dim.ndime; npnod=dim.npnod; nunkn = dim.nunkn; nstre = dim.nstre;

ptype = problembsc.problemtype;
ftype = problembsc.phisical_type;

switch ftype
        case {'ELASTIC'}
               switch element.material.homogenous_material
                    case 'YES'
                        switch element.material.type
                            case 'HOOKE_LAW'
                                mu = element.mu_func(gamma);
                                kappa = element.kappa_func(gamma);
                                
%                                 epoiss = (3*kappa - 2*mu)./(2*(3*kappa+mu));
                                epoiss = (kappa - mu)./(kappa + mu);
                                eyoung = 4*kappa.*mu./(kappa + mu);
%                                 eyoung = (9*kappa.*mu)./(3*kappa+mu);
                                
                                
                                %epoiss = ones(1,nelem)*element.material.poiss;
                                %eyoung = ones(1,nelem)*element.material.young;
                                
                                [Ce] = hooke_law_celas(nstre,nelem,element,ptype,eyoung,epoiss,igaus);
%                                 Ce_plus = Ce;
%                                 Ce_minus = element.material.opt_epsi*Ce;
%                                 
%                                 Ce(1,1,:) = avarage_element(squeeze(Ce_plus(1,1,:)),squeeze(Ce_minus(1,1,:)),chi');
%                                 Ce(1,2,:) = avarage_element(squeeze(Ce_plus(1,2,:)),squeeze(Ce_minus(1,2,:)),chi');
%                                 Ce(1,3,:) = avarage_element(squeeze(Ce_plus(1,3,:)),squeeze(Ce_minus(1,3,:)),chi');
%                                 Ce(2,1,:) = avarage_element(squeeze(Ce_plus(2,1,:)),squeeze(Ce_minus(2,1,:)),chi');
%                                 Ce(2,2,:) = avarage_element(squeeze(Ce_plus(2,2,:)),squeeze(Ce_minus(2,2,:)),chi');
%                                 Ce(2,3,:) = avarage_element(squeeze(Ce_plus(2,3,:)),squeeze(Ce_minus(2,3,:)),chi');
%                                 Ce(3,2,:) = avarage_element(squeeze(Ce_plus(3,2,:)),squeeze(Ce_minus(3,2,:)),chi');
%                                 Ce(3,3,:) = avarage_element(squeeze(Ce_plus(3,3,:)),squeeze(Ce_minus(3,3,:)),chi');
                                
                                
                                
                        end
                        
                    case 'VADEMECUM'
                        Ch = squeeze(element.Ch(igaus,:,:));
                        
                        Ce_plus(1,1,:) = Ch(1,:);
                        Ce_plus(1,2,:) = Ch(2,:);
                        Ce_plus(1,3,:) = Ch(3,:);
                        Ce_plus(2,1,:) = Ch(2,:);
                        Ce_plus(2,2,:) = Ch(4,:);
                        Ce_plus(2,3,:) = Ch(5,:);
                        Ce_plus(3,2,:) = Ch(5,:);
                        Ce_plus(3,3,:) = Ch(6,:);
                        
                        Ce_minus = element.material.opt_epsi*Ce_plus;
                        
                        Ce(1,1,:) = avarage_element(squeeze(Ce_plus(1,1,:)),squeeze(Ce_minus(1,1,:)),chi');
                        Ce(1,2,:) = avarage_element(squeeze(Ce_plus(1,2,:)),squeeze(Ce_minus(1,2,:)),chi');
                        Ce(1,3,:) = avarage_element(squeeze(Ce_plus(1,3,:)),squeeze(Ce_minus(1,3,:)),chi');
                        Ce(2,1,:) = avarage_element(squeeze(Ce_plus(2,1,:)),squeeze(Ce_minus(2,1,:)),chi');
                        Ce(2,2,:) = avarage_element(squeeze(Ce_plus(2,2,:)),squeeze(Ce_minus(2,2,:)),chi');
                        Ce(2,3,:) = avarage_element(squeeze(Ce_plus(2,3,:)),squeeze(Ce_minus(2,3,:)),chi');
                        Ce(3,2,:) = avarage_element(squeeze(Ce_plus(3,2,:)),squeeze(Ce_minus(3,2,:)),chi');
                        Ce(3,3,:) = avarage_element(squeeze(Ce_plus(3,3,:)),squeeze(Ce_minus(3,3,:)),chi');
                        
                        
                        %epoiss = element.material.epoiss_micro; 
                end    
                    
            
        case {'THERMAL'}
            switch element.material.type
                case 'FOURIER_LAW'
                    kxx = ones(1,nelem).*element.material.k11;
                    kxy = ones(1,nelem).*element.material.k12;
                    kyy = ones(1,nelem).*element.material.k22;
                    [Ce_plus] = fourier_law_elas(chi,kxx,kxy,kyy);
                    Ce_minus = element.material.opt_epsi*Ce_plus;
                    
                    Ce(1,1,:) = avarage_element(squeeze(Ce_plus(1,1,:)),squeeze(Ce_minus(1,1,:)),chi');
                    Ce(1,2,:) = avarage_element(squeeze(Ce_plus(1,2,:)),squeeze(Ce_minus(1,2,:)),chi');
                    Ce(2,1,:) = avarage_element(squeeze(Ce_plus(2,1,:)),squeeze(Ce_minus(2,1,:)),chi');
                    Ce(2,2,:) = avarage_element(squeeze(Ce_plus(2,2,:)),squeeze(Ce_minus(2,2,:)),chi');
 
                    
                    
                    
            end
end
end
    
    
 