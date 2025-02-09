function U = unionHull(varargin)
%unionHull - computes the convex hull of the union of other bounded convex polyhedra. 
%Any number of polyhedra can be input to the union operation. They can
%be expressed using either vertices or linear in/equalities, according to
%the input scheme,
%
%
%   I = intersectionHull('vert', V1, 'lcon', A2,b2,...,'qlcon', x3,A3,b3,...)
%
%The arguments specifying different polyhedra are separated using labels
%'vert', 'lcon', or 'qlcon'. The label 'vert' signifies that an input polyhedron will
%be expressed using vertices and is to be followed by any string of input arguments 
%accepted by vert2lcon(). The label 'lcon' signifies that an input polyhedron will
%be expressed using linear constraints and is to be followed by any string of 
%input arguments accepted by lcon2vert(). The label 'qlcon' signifies that an input 
%polyhedron will be expressed using linear constraints and a known interior point and
%is to be followed by  input arguments accepted by qlcon2vert().
%
%The output, U, is a struct containing fields
%
%   U.vert: A matrix whose rows are the vertices of the output polyhedron.
%   U.lcon: The quadruplet of linear constraint data {A,b,Aeq,beq}
%           describing the output polyhedron.
%
%EXAMPLE 1: The unit simplex in 3D is formed as the convex hull of the union 
%of a cube and a triangle, both expressed in terms of their vertices.
% 
%     V1=(dec2bin(0:2^3-1,3)-'0')/3.0001;    %vertices of 3D cube
%     V2=eye(3);  %vertices of triangle
% 
%     U=unionHull('vert',V1,'vert',V2);   %compute the union
%
%One can verify that the vertices of the result are those of the unit
%simplex,
%
%       >> U.vert
% 
%             ans =
% 
%                  0     0     0
%                  0     0     1
%                  0     1     0
%                  1     0     0
%
%We also demonstrate the use of separateBounds() to convert the result to 
%Optimization Toolbox representation. The tolerance parameter 0.99 is used
%to deal with residual floating point noise.
%
%         >> [A,b,Aeq,beq,lb,ub]=separateBounds(U.lcon{:},.99)
% 
%             A =
% 
%                 0.5774    0.5774    0.5774
% 
% 
%             b =
% 
%                 0.5774
%
% 
%             Aeq =
% 
%                  []
% 
% 
%             beq =
% 
%                Empty matrix: 0-by-1
% 
% 
%             lb =
% 
%                  0
%                  0
%                  0
% 
% 
%             ub =
% 
%                Inf
%                Inf
%                Inf   
%   
% 
%   
%EXAMPLE 2: The unit cube is formed as the union of a simplex, a square,
%and a single 3D point. The simplex and square are both represented using
%equalites and inequalities. For convenience, we use the addBounds() utility
%to compose the in/equality matrices. The representation of the square illustrates
%the use of the 'qlcon' label, when a known interior point is available.
% 
%   
%     [A1,b1,Aeq1,beq1] = addBounds([1 1 1],1,[],[],[0;0;0]); %unit simplex
% 
%     [A2,b2,Aeq2,beq2] = addBounds([],[],[],[],[1;0;0],[1;1;1]); %a square face of the unit cube
% 
%                    x2 = [1,0.5,0.5];%a point in the relative interior of the square  
% 
%     V3 = [0,1,1]; %one vertex of the unit cube
% 
% 
%     U=unionHull('lcon',A1,b1,Aeq1,beq1,...         %compute the union
%                 'qlcon',x2,A2,b2,Aeq2,beq2,...
%                  'vert',V3);
%
%One can verify that the result contains the 8 vertices of the unit cube,
%
%       >> U.vert
% 
%         ans =
% 
%                  0   -0.0000   -0.0000
%             0.0000   -0.0000    1.0000
%            -0.0000    1.0000   -0.0000
%                  0    1.0000    1.0000
%             1.0000    0.0000    0.0000
%             1.0000    0.0000    1.0000
%             1.0000    1.0000    0.0000
%             1.0000    1.0000    1.0000
%
%As before, we can use separateBounds() to convert to Optimization Toolbox
%form,
% 
%     >>  [A,b,Aeq,beq,lb,ub]=separateBounds(U.lcon{:},.99)  
%
%which results in A=b=Aeq=beq=[] and
%
%    lb =
% 
%          0
%          0
%          0
% 
% 
%     ub =
% 
%         1.0000
%         1.0000
%         1.0000


%%%%begin parsing

if isnumeric(varargin{1})
   TOLcell=varargin(1);
   varargin(1)=[];
else
    TOLcell={};
end

N=length(varargin);
idxType = [find(cellfun(@ischar,varargin)),N+1];

L=length(idxType)-1;
S(L).type=[];
S(L).args={};


for i=1:L
   
    j=idxType(i);
    k=idxType(i+1);
    S(i).type=varargin{j};
    S(i).args=varargin(j+1:k-1);
    
    if isempty(S(i).args)
     error 'Syntax error - arguments missing' 
    end
    
     switch S(i).type
        
        case 'lcon'
            
             S(i).V = lcon2vert(S(i).args{:});
            
        case 'qlcon'
            
            S(i).V = qlcon2vert(S(i).args{:});
            
        case 'vert'
            
            S(i).V = S(i).args{1};
            
        otherwise
            
            error(['Unrecognized representation label of polyhedron ' num2str(i)]);
        
        
      end
    

 end



%%%%end parsing


    V=vertcat(S.V);
   
    x0=mean(V,1);
    
   [U.lcon{1:4}]=vert2lcon(V,TOLcell{:});
   
   U.vert=qlcon2vert(x0,U.lcon{:},TOLcell{:});
  

   
