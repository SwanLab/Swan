function [Ce] = hooke_law_celas(nstre,nelem,element,ptype,eyoung,epoiss,igaus)

Ce = zeros(nstre,nstre,nelem);
switch ptype
    case '2D'
        switch element.material.subtype
            case 'PLANESTRAIN'
                %c1 = chi.*eyoung./((1+epoiss).*(1-2*epoiss));
                c1 = eyoung./((1+epoiss).*(1-2*epoiss));
                Ce(1,1,:) = c1.*(1-epoiss);
                Ce(1,2,:) = c1.*epoiss;
                Ce(2,1,:) = c1.*epoiss;
                Ce(2,2,:) = c1.*(1-epoiss);
                Ce(3,3,:) = c1*0.5.*(1-2*epoiss);
                
            case 'PLANESTRES'
                
                switch element.material.homogenous_material
                    case 'YES'
                        %c1 = eyoung./(1-epoiss.^2);
                        % stres(1,:) = c1.*(strain(1,:) + epoiss.*strain(2,:));   % sigma11
                        %  stres(2,:) = c1.*(epoiss.*strain(1,:) + strain(2,:));   % sigma22
                        % stres(3,:) = c1.*(0.5*(1-epoiss).*strain(3,:));         % sigma12
                        % strain(4,:) = -epoiss.*(stres(1,:)+stres(2,:))./eyoung;  % strain33
                        %c1 = chi.*eyoung./(1-epoiss.^2);
                        
                        c1 = eyoung./(1-epoiss.^2);                        
                        Ce(1,1,:) = c1;
                        Ce(1,2,:) = c1.*epoiss;
                        Ce(2,1,:) = c1.*epoiss;
                        Ce(2,2,:) = c1;
                        Ce(3,3,:) = c1*0.5.*(1-epoiss);
                        
%                     case 'VADEMECUM'
%                         Ch = bsxfun(@times,squeeze(element.Ch(igaus,:,:)),chi);
%                         Ce(1,1,:) = Ch(1,:);
%                         Ce(1,2,:) = Ch(2,:);
%                         Ce(1,3,:) = Ch(3,:);
%                         Ce(2,1,:) = Ch(2,:);
%                         Ce(2,2,:) = Ch(4,:);
%                         Ce(2,3,:) = Ch(5,:);
%                         Ce(3,2,:) = Ch(5,:);
%                         Ce(3,3,:) = Ch(6,:);
                        
                end
     
        end
    case '3D'
        c0 = eyoung./((1+epoiss).*(1-2*epoiss));
        %c0 = chi.*eyoung./((1+epoiss).*(1-2*epoiss));
        c11 = c0.*(1-epoiss);
        c12 = c0.*epoiss;
        emu = eyoung./(2*(1+epoiss));
        %emu = chi.*eyoung./(2*(1+epoiss));
        Ce(1,1,:) = c11;
        Ce(1,2,:) = c12;
        Ce(1,3,:) = c12;
        
        Ce(2,1,:) = c12;
        Ce(2,2,:) = c11;
        Ce(2,3,:) = c12;
        
        Ce(3,1,:) = c12;
        Ce(3,2,:) = c12;
        Ce(3,3,:) = c11;
        
        Ce(4,4,:) = emu;
        Ce(5,5,:) = emu;
        Ce(6,6,:) = emu;
end

end

