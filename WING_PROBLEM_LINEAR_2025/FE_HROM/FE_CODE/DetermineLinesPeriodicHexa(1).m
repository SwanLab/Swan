function NODESln = DetermineLinesPeriodic(NODESpl)


%

if nargin == 0
    load('tmp.mat')
end




NODESln = cell(size(NODESpl,1),size(NODESpl,2),3,2) ;
for idimA = 1:size(NODESpl,1)
     for jdimA = 1:size(NODESpl,2)    
         NODA = NODESpl{idimA,jdimA} ; 
         for idimB =idimA+1:size(NODESpl,1)
             for jdimB = 1:size(NODESpl,2)    
                 NODB = NODESpl{idimB,jdimB} ; 
                     NODESln{idimA,jdimA,idimB,jdimB} = intersect(NODA,NODB);
                 
%                    
                 
             end
         end
         
   
    end
end

%
% for iplane = 1:length(PLANESf)
%     iname = PLANESf{iplane} ;
%     for jplane = 1:length(PLANESf)
%         jname = PLANESf{jplane} ;
%         ijname = [iname,'_',jname] ;
%         iijj = intersect(NODESpl.(iname),NODESpl.(jname));
%         if  jplane>iplane && ~isempty(iijj) && length(iijj) ~= length(NODESpl.(iname))
%
%         NODESln = setfield(NODESln,ijname,iijj) ;
%
%         end
%     end
% end