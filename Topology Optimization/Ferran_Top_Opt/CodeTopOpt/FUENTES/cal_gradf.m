function [gradF] = cal_gradf(problemtype,coordn,coorda,cartd,nelem,gradF,igaus,gpxy)
gi=2/3; gj=1/6; 
posgp=[gi gj gj;
       gj gi gj];

gradF = zeros(3,3,nelem);
for i=1:2
    for j=1:2
        for a=1:3
            gradF(i,j,:)=gradF(i,j,:)+coorda(a,i,:).*cartd(j,a,:);
        end
    end
end

switch problemtype
    case 'PLANESTRAIN'
        gradF(3,3,:) = ones(nelem,1);
    case 'AXI'
        %calculo del desplazamiento en el pto de gauss
        exisp = posgp(1,igaus);
        etasp = posgp(2,igaus);
        shape = shapef(exisp,etasp);
        rdisp = zeros(2,nelem);
%        keyboard
        for idime=1:2
            for inode=1:3
                rdisp(idime,:) = rdisp(idime,:)+...
                squeeze((coorda(inode,idime,:)-coordn(inode,idime,:)))'.*shape(inode);
            end
        end
        gradF(3,3,:)=rdisp(1,:)./squeeze(gpxy(igaus,1,:))'+1;
        
%         radio=(coordn(1,1,:)+coordn(2,1,:)+coordn(3,1,:))/3;
%         rdisp=((coorda(1,1,:)-coordn(1,1,:))+(coorda(2,1,:)-coordn(2,1,:))+...
%             (coorda(3,1,:)-coordn(3,1,:)))/3;
%         gradF(3,3,:)=rdisp./radio+1;

end
        
end

