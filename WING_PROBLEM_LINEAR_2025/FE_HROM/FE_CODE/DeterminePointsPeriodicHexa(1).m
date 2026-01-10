function NODESpt = DeterminePointsPeriodicHexa(NODESpl)

if nargin == 0
    load('tmp.mat')
end


NODESpt = [] ;
for idimA = 1:size(NODESpl,1)
     for jdimA = 1:size(NODESpl,2)
         NODA = NODESpl{idimA,jdimA} ; 
         for idimB =idimA+1:size(NODESpl,1) 
             for jdimB = 1:size(NODESpl,2) 
                 NODB = NODESpl{idimB,jdimB} ; 
                PPP =   intersect(NODA,NODB);      
                 for  idimC = idimB+1:size(NODESpl,1)
                     for  jdimC = 1:size(NODESpl,2)
                                               NODC = NODESpl{idimC,jdimC} ; 

                    NODESpt = [ NODESpt intersect(PPP,NODC)];       
                     end
                 end
             end
         end
         
   
    end
end
NODESpt = NODESpt' ; 
% 
% 
% PLANESf = fieldnames(NODESpl);
% 
% NODESpnt = [] ;
% for iplane = 1:length(PLANESf)
%     iname = PLANESf{iplane} ;
%     for jplane = 1:length(PLANESf)
%         jname = PLANESf{jplane} ;
%         ijname = [iname,'_',jname] ;
%         iijj = intersect(NODESpl.(iname),NODESpl.(jname));
%         if  jplane>iplane && ~isempty(iijj) && length(iijj) ~= length(NODESpl.(iname))
%             for kplane =    1:length(PLANESf)
%                 %NODESln = setfield(NODESln,ijname,iijj) ;
%                 kname = PLANESf{kplane} ;
%                 ijkname = [iname,'_',jname,'_',kname] ;
%                 iijjkk = intersect(iijj,NODESpl.(kname));
%                 if  kplane>jplane &&  length(iijjkk)==1
%                     
%                     NODESpnt = setfield(NODESpnt,ijkname,iijjkk) ;
%                     
%                 end
%                 
%                 
%             end
%             
%         end
%     end
% end