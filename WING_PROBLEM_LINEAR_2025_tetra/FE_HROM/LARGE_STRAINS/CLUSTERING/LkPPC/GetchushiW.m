function   W=GetchushiW(cX,k,c)
 [~,n6]=size(cX);
   for i=1:k
           tA=cX((cX(:,1)==i),2:n6);
           tB=cX((cX(:,1)~=i),2:n6);
          W(i,:)=GepOneSide(tA,tB,c);
  end         
end
