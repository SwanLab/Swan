function I = intersectionHull(varargin)
%intersectionHull - computes a bounded convex polyhedron resulting from the
%intersection of other convex polyhedra. The other polyhedra need not individually
%be bounded. Only the final intersection resulting from them must be bounded. 
%Any number of polyhedra can be input to the intersection operation. They can
%be expressed using either vertices or linear in/equalities, according to
%the input scheme,
%
%
%   I = intersectionHull('vert', V1, 'lcon', A2,b2, 'lcon', A3,b3,Aeq3,,beq3,...)
%
%The arguments specifying different polyhedra are separated using labels
%'vert' and 'lcon'. The label 'vert' signifies that an input polyhedron will
%be expressed using vertices and is to be followed by any string of input arguments 
%accepted by vert2lcon(). The label 'lcon' signifies that an input polyhedron will
%be expressed using linear constraints and is to be followed by any string of 
%input arguments accepted by lcon2vert().
%
%The output, I, is a struct containing fields
%
%   I.vert: A matrix whose rows are the vertices of the polyhedron formed from 
%           the intersection.
%   I.lcon: The quadruplet of linear constraint data {A,b,Aeq,beq}
%           describing the polyhedral intersection.
%
%EXAMPLE 1: This example computes the intersection of a unit square and an 
%oblique 2D line segment, both expressed in terms of their vertices.
% 
%     V1=dec2bin(0:2^2-1,2)-'0';   %vertices of unit square
% 
%     V2=[1,1;0,-1];               %vertices of 2D line segment
% 
%     I=intersectionHull('vert',V1,'vert',V2);  %compute intersection
%
%The intersection is another line segment with vertices
%
%        >> I.vert  
% 
%         ans =
% 
%             0.5000         0
%             1.0000    1.0000 
% 
%EXAMPLE 2: This example computes the intersection of a unit cube, expressed in
%terms of its vertices, and an infinite oblique 3D line, expressed in terms of linear equalities. 
%Note that the line is an unbounded polyhedron. This is okay, since we know in advance 
%that the final polyhedron formed from the intersection is bounded.
% 
%     V=dec2bin(0:2^3-1,3)-'0';    %vertices of unit cube
% 
%     Aeq=[1 -1 0; 0 1 -1]; beq=[0;0]; %oblique line in 3D
% 
%     I=intersectionHull('vert',V,'lcon',[],[],Aeq,beq);  %compute intersection
% 
%Once again, the intersection is a line segment. Its vertices are
%
%      >> I.vert   %vertices of line segment of intersection
%
%         ans =
% 
%             0.0000    0.0000    0.0000
%             1.0000    1.0000    1.0000


%%%%begin parsing

if isnumeric(varargin{1})
   TOL=varargin{1};
   varargin(1)=[];
else
    TOL=[];
end

N=length(varargin);
idxType = [find(cellfun(@ischar,varargin)),N+1];

L=length(idxType)-1;
S(L).type=[];
S(L).args={};
S(L).A=[];
S(L).b=[];
S(L).Aeq=[];
S(L).beq=[];

for i=1:L
   
    j=idxType(i);
    k=idxType(i+1);
    S(i).type=varargin{j};
    S(i).args=varargin(j+1:k-1);
    
    if isempty(S(i).args)
     error 'Syntax error - arguments missing' 
    end
    
    lcon=cell(1,4);
    
      switch S(i).type
  
          case 'vert'
              
              [lcon{1:4}] = vert2lcon(S(i).args{:}); 
              
          case 'lcon'
              
               lcon(1:k-j-1) = S(i).args;
              
          case 'qlcon' %deliberately undocumented - no point in using this
              
               lcon(1:k-j-1) = S(i).args(2:end);
               
          otherwise
              
            error(['Unrecognized representation label of polyhedron ' num2str(i)]);
          
      end

    
    [S(i).A, S(i).b, S(i).Aeq, S(i).beq] = deal(lcon{:});
    
end


%%%%end parsing


   A=vertcat(S.A);
   b=vertcat(S.b); 
   Aeq=vertcat(S.Aeq);
   beq=vertcat(S.beq); 
   
   [V,nr,nre]=lcon2vert(A,b,Aeq,beq,TOL);
   
   I.vert=V;
   I.lcon={A(nr,:),b(nr,:), Aeq(nre,:),beq(nre,:)};

   
