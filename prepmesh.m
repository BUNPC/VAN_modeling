function [ii,jj,connnum,barea,type1node,freenode]=prepmesh(node,elem,boundary)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   input: node, elem, boundary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% nn: number of nodes, ne: number of elements
nn=length(node);
ne=length(elem);

% idx&bb: all surf. elem with tag 0 (vessel surface)
% type1node: indices of Type-1 BC nodes (vessel surface)
% freenode:  indices of the remaining nodes

idx=find(boundary(:,end)==0);
bb=boundary(idx,1:3);
allnodes=zeros(nn,1);
allnodes(bb(:))=1;

type1node=find(allnodes==1);
freenode=find(allnodes==0);

barea=elemvolume(node,boundary(idx,1:3));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Calculate element volume & the indices of non-zero
%   elements of LHS matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(1,'analyzing mesh connections ...\n');
tic;
vol=elemvolume(node,elem);
[ii,jj,connnum]=femnz(elem,nn);
toc

