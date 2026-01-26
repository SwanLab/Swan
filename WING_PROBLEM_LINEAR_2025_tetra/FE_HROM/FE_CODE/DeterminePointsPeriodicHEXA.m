function[NODESpt LINES]= DeterminePointsPeriodicHEXA(NODESln)

if nargin == 0
    load('tmp.mat')
end


NODESpt = zeros(12,1) ;  
LINES = zeros(12,2) ; 

LINES(1,:) = [1 2] ; 
LINES(2,:) = [2 3] ; 
LINES(3,:) = [3 4] ; 
LINES(4,:) = [4 5] ; 
LINES(5,:) = [5 6] ; 
LINES(6,:) = [6 1] ; 
LINES(7,:) = [7 8] ; 
LINES(8,:) = [8 9] ; 
LINES(9,:) = [9 10] ; 
LINES(10,:) = [10 11] ; 
LINES(11,:) = [11 12] ; 
LINES(12,:) = [12 7] ; 

for i = 1:size(LINES,1) 
    
    NODESpt(i) = intersect(NODESln{LINES(i,1)},NODESln{LINES(i,2)}) ; 
    
end



 
 



% 
% NODESpt = [] ;
% for idimA = 1:size(NODESpl,1)
%     for jdimA = 1:size(NODESpl,2)
%         NODA = NODESpl{idimA,jdimA} ;
%         for idimB =idimA+1:size(NODESpl,1)
%             for jdimB = 1:size(NODESpl,2)
%                 NODB = NODESpl{idimB,jdimB} ;
%                 PPP =   intersect(NODA,NODB);
%                 for  idimC = idimB+1:size(NODESpl,1)
%                     for  jdimC = 1:size(NODESpl,2)
%                         NODC = NODESpl{idimC,jdimC} ;
%                         
%                         NODESpt = [ NODESpt intersect(PPP,NODC)];
%                     end
%                 end
%             end
%         end
%         
%         
%     end
% end
% NODESpt = NODESpt' ;
% %
% %
% % PLANESf = fieldnames(NODESpl);
% %
% % NODESpnt = [] ;
% % for iplane = 1:length(PLANESf)
% %     iname = PLANESf{iplane} ;
% %     for jplane = 1:length(PLANESf)
% %         jname = PLANESf{jplane} ;
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