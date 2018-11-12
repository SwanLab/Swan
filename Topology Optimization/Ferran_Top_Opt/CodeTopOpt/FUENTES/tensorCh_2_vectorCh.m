function Ch = tensorCh_2_vectorCh(Ch_tens)

Ch(:,1) = Ch_tens(1,1,:);  
Ch(:,2) = Ch_tens(1,2,:);  
Ch(:,3) = Ch_tens(1,3,:);  
Ch(:,2) = Ch_tens(2,1,:);  
Ch(:,4) = Ch_tens(2,2,:);  
Ch(:,5) = Ch_tens(2,3,:);  
Ch(:,3) = Ch_tens(3,1,:);  
Ch(:,5) = Ch_tens(3,2,:);  
Ch(:,6) = Ch_tens(3,3,:);  


end
