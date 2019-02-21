function VAN_setup_mesh(mousefolder1,savefolder,roi)

%%-------------------------------------------------------------------------
% Inputs:
%       roi: Structure with the ROI dimensions in um
%            roi.xx = [xmin xmax] (um)
%            roi.yy = [ymin ymax] (um)
%            roi.zz = [zmin zmax] (um)
%
%       opt: Parameters for the vol2mesh function. Copied from 'vol2surf' help:
%            opt = a float number>1: max radius of the Delaunay sphere(element size) 
%            opt.radbound: same as above, max radius of the Delaunay sphere
%            opt.distbound: maximum deviation from the specified isosurfaces, default is 1
%            opt(1,2,...).radbound: same as above, for each levelset
%
%--------------------------------------------------------------------------


%%--------------------------------------------------------------------------
% We will generate a mesh but only on a subportion of the stack to avoid
% regions with plugged capillaries and therefore no oxygen reaching that
% point. 

% For that we will need a sub-graph containing only segments entirely
% conainted inside of the ROI (otherwise we would be creating extra BCs).
% BUT THE FLOW MUST CONSERVED OTHERWISE PO2 WILL BE PILLING UP ! So we need
% to consider an extra output term each time we remove a segment to make
% sure the oxygen is not pilling up. What about when the removed segment
% was inflowing to that node? How do we set the O2 amount inflowing? This
% is a major problem and we need another solution. 

% Another way is to keep the full graph for the O2 advection on the graph.
% Then, we run the O2 diffusion on the ROI volume. The portion of the graph
% outside of the diffusion volume will be impermeable for now. We can
% modify that later to give them an average permeability per unit length
% since there is not tissue arround them and that we don't know the
% O2 gradient from inside to outside. 
% Let's try that.
%--------------------------------------------------------------------------


%%--------------------------------------------------------------------------
% created 1/27/2013 by L. Gagnon
% edited 08/31/2016 by J. Selb to make it a function
%--------------------------------------------------------------------------



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define ROI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([mousefolder1,'graph.seed'], '-mat')

if isfield(im2,'Vm_1um')
opt.radbound=30;
opt.distbound=1; %default is 1%% in the file it says 0.8
ROIx = roi.xx;    %um
ROIy = roi.yy;    %um
ROIz = roi.zz;    %um

xidx1=max(round(ROIx(1)),1);
xidx2=min(round(ROIx(2)),size(im2.Vm_1um,1));
yidx1=max(round(ROIy(1)),1);
yidx2=min(round(ROIy(2)),size(im2.Vm_1um,2));
zidx1=max(round(ROIz(1)),1);
zidx2=min(round(ROIz(2)),size(im2.Vm_1um,3));

%isotropic 1um mask
img = im2.Vm_1um(xidx1:xidx2,yidx1:yidx2,zidx1:zidx2);
else
  opt.radbound=30;
opt.distbound=1; %default is 1%% in the file it says 0.8
ROIx = roi.xx;    %um
ROIy = roi.yy;    %um
ROIz = roi.zz;    %um


xidx1=max(round(ROIx(1)),1);
xidx2=min(round(ROIx(2)),size(im2.Vm,1));
yidx1=max(round(ROIy(1)),1);
yidx2=min(round(ROIy(2)),size(im2.Vm,2));
zidx1=max(round(ROIz(1)),1);
zidx2=min(round(ROIz(2)),size(im2.Vm,3));




% xidx1=max(round(ROIx(1)),1);                          %xidx1=round(xidx1/im2.Hvox(1));
% xidx2=min(round(ROIx(2)),floor(size(im2.Vm,1)*im2.Hvox(1)));         %xidx2=round(xidx2/im2.Hvox(1));
% yidx1=max(round(ROIy(1)),1);                           %yidx1=round(yidx1/im2.Hvox(2));
% yidx2=min(round(ROIy(2)),floor(size(im2.Vm,2)*im2.Hvox(2)));          %yidx2=round(yidx2/im2.Hvox(2));
% zidx1=max(round(ROIz(1)),1);                           %zidx1=round(zidx1/im2.Hvox(3));
% zidx2=min(round(ROIz(2)),floor(size(im2.Vm,3)*im2.Hvox(3)));          %zidx2=round(zidx2/im2.Hvox(3));

%isotropic 1um mask
img = im2.Vm(xidx1:xidx2,yidx1:yidx2,zidx1:zidx2);  
im2.Vm_1um=im2.Vm;
    
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mesh the reduce stack
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[node,elem,boundary] = v2m( (img==0), 0.5, opt, 100, 'cgalmesh'); % 200 000 nodes

clear img; %to save memory

% transpose x and y to be consistent with matlab
node = node(:,[2 1 3]); 



% Now we need to set up boundary to distinguish
% bounding box from vessel wall node

% get bounding box corners
minval = min(node);
maxval = max(node);

% compute the centroid of all surface triangles
bc = meshcentroid(node,boundary(:,1:3));
boundary(:,end) = 0;
tol = (maxval(1)-minval(1))/100;
boundary( find(bc(:,1)<minval(1)+tol | bc(:,1)>maxval(1)-tol | ...
    bc(:,2)<minval(2)+tol | bc(:,2)>maxval(2)-tol | ...
    bc(:,3)<minval(3)+tol | bc(:,3)>maxval(3)-tol ),end) = 1;
[iiPrepMesh,jjPrepMesh,connnum,barea,vesWallnode,freenode]=prepmesh(node(:,1:3),elem(:,1:4),boundary);

% visualize the vessels
if 1
%    foo=ismember(elem(:,1:3),vesWallnode);
%    lst = find(sum(foo,2)==3);
    figure;
    trisurf( boundary(find(boundary(:,end)==0),1:3), node(:,1), node(:,2), -node(:,3), 'linestyle','none','facecolor','flat');
    
end    
  
nn=size(node,1);
nng = size(im2.nodePos,1);
nodearea=zeros(nn,1);
idx=find(boundary(:,end)==0);
for i=1:length(idx)
    nodearea(boundary(idx(i),1:3))=nodearea(boundary(idx(i),1:3))+barea(i);
end
nodearea=nodearea/3.0;
vol=elemvolume(node(:,1:3),elem(:,1:4));
nodevol=zeros(nn+nng,1);
for i=1:size(elem,1)
    nodevol(elem(i,1:4))=nodevol(elem(i,1:4))+vol(i);
end
nodevol=nodevol/4.0;

%re-adjust the offset (x and y inverted)
node(:,1)=node(:,1)+(yidx1-1);
node(:,2)=node(:,2)+(xidx1-1);
node(:,3)=node(:,3)+(zidx1-1);

  

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate the vol of each node given vessel Mask
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
volg = getNodeVol( im2.Vm_1um, im2.nodePos_um );
volg = max(volg,1); % This needs a better solution. I suspect I am getting 
                    % elements of volg=0 because the node is not in the Vm
                    % mask

                           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make advection matrix 
% a positive flow indicate flow from node in col 1
% to node in col2
% a negative flow goes in the other direction
%
% \Delta{n} / dt    = vel \dot \grad{n} 
% vol \Delta{c} /dt = vel \dot \grad(vol c}
%                   = (vel area) len (\Delta{c} / len)
%                   = flow \Delta{c}
% \Delta{c} / dt    = flow \ Delta{c} / vol
%
% ( F_{i-1} C_{i-1}^k - F_i C_i^k ) dt = V_i ( C_i^{k+1} - C_i^k )
%
% INPUTS
%   volg - volume associated with each node. Calculated by
%          volg = getNodeVol( im2.Vm, im2.nodePos );
%   im2.nodeEdges
%   im2.edgeFlow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
neg = size(im2.nodeEdges,1);
nng = size(im2.nodePos,1);

Ladv = zeros(nng,nng);
nodeIn = zeros(nng,1);
Fnode = zeros(nng,1);

LInidx = zeros(nng,4);            % Ladv indices for in flow
LInFeidx = (neg+1)*ones(nng,4);          % edge indices for in flow
LOutidx = zeros(nng,1);           % Ladv indices for out flow
LOutFeidx = (neg+1)*ones(nng,4);  % edge indices for out flow

Fedges = im2.edgeFlow;
Fedges(neg+1) = 0;                %   use neg+1 index for 0 flow
im2.edgeFlow = Fedges(1:neg+1);

eg = im2.nodeEdges;

for ii=1:nng
    [ir,ic] = find( eg==ii );
    Fnode(ii) = mean(Fedges(ir));
    for jj=1:length(ir)
        n1 = eg(ir(jj),1);
        n2 = eg(ir(jj),2);
        eIdx = ir(jj);
        if (n2==ii && Fedges(eIdx)>=0) 
            Ladv(ii,n1) = abs(Fedges(eIdx)) / volg(ii);
            LInidx(ii,jj) = ii + (n1-1)*nng;
            LInFeidx(ii,jj) = eIdx;
        elseif (n1==ii && Fedges(eIdx)<0)
            Ladv(ii,n2) = abs(Fedges(eIdx)) / volg(ii);
            LInidx(ii,jj) = ii + (n2-1)*nng;
            LInFeidx(ii,jj) = eIdx;
        elseif (n2==ii && Fedges(eIdx)<0) || (n1==ii && Fedges(eIdx)>=0)
            Ladv(ii,ii) = Ladv(ii,ii) -1*abs(Fedges(eIdx)) / volg(ii);
            LOutidx(ii) = ii + (ii-1)*nng;
            LOutFeidx(ii,jj) = eIdx;
        end
    end
    if length(ir)==1 && ((n2==ii && Fedges(eIdx)>=0) || (n1==ii && Fedges(eIdx)<0))
        Ladv(ii,ii) =  -1*abs(Fedges(eIdx)) / volg(ii);
        nodeIn(ii) = -1;
        if n1==ii;
            nn = n2;
        else
            nn = n1;
        end
        LInidx(ii,1) = ii + (nn-1)*nng;
        LInFeidx(ii,1) = eIdx;
        LOutidx(ii) = ii + (ii-1)*nng;
        LOutFeidx(ii,1) = eIdx;
    elseif length(ir)==1 && ((n2==ii && Fedges(eIdx)<0) || (n1==ii && Fedges(eIdx)>=0))
        Ladv(ii,ii) = 0;
        nodeIn(ii) = 1;
        LInidx(ii,1) = ii + (ii-1)*nng;
        LInFeidx(ii,1) = neg+1;
        LOutidx(ii) = ii + (ii-1)*nng;
        LOutFeidx(ii,1) = neg+1;
    end
end
Ladv = sparse(Ladv);



                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% map each vessel wall node to a graph node
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Note: Because of the ROI, each vessel wall node will be mapped to a graph
%node but some graph node will be left without a vessel wall node. This is
%OK as advection is run on graph node while diffusion is run on vessel wall
%nodes independently.

ng = im2.nodePos_um;
gfMap = zeros(nn,1);
for idxN = 1:length(vesWallnode)
    r2 = ( (ng(:,1)-node(vesWallnode(idxN),1)).^2 + ...
          (ng(:,2)-node(vesWallnode(idxN),2)).^2 + ...
          (ng(:,3)-node(vesWallnode(idxN),3)).^2 );
    [foo,idx] = min(r2);
    gfMap(vesWallnode(idxN)) = idx(1);
    
%     if sqrt(foo) > 20
%         display(sprintf('elem: %d and node: %d were %2.1f voxels appart',vesWallnode(idxN),idx(1),sqrt(foo)))
%     end
        
end



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make flux matrix (now coded using sparse matrix LG)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Note: This needs to be changed because of the ROI. What do I do with the
%graph nodes outside of the mesh volume? I will just skip them for now and
%see what I get.

% Lflux = zeros(nng,length(vesWallnode));
Lflux = sparse([],[],[],nng,length(vesWallnode)); %defined sparse to save memory
for idxG = 1:nng
    lst = find(gfMap(vesWallnode)==idxG);
    if ~isempty(lst)
        Lflux(idxG,lst) = 1;
    end
end
% Lflux = sparse(Lflux);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save in data structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Store mesh info
im2.Mesh.node = node(:,1:3);
im2.Mesh.elem = elem(:,1:4);
im2.Mesh.boundary = boundary;
im2.Mesh.iiPrepMesh = iiPrepMesh;
im2.Mesh.jjPrepMesh = jjPrepMesh;
im2.Mesh.connnum = connnum;
im2.Mesh.vol = vol;
im2.Mesh.nodearea = nodearea;
im2.Mesh.nodevol = nodevol;
im2.Mesh.vesWallnode = vesWallnode;

% store advection info
im2.Adv.Ladv = Ladv;
im2.Adv.LInidx = LInidx;
im2.Adv.LInFeidx = LInFeidx;
im2.Adv.LOutidx = LOutidx;
im2.Adv.LOutFeidx = LOutFeidx;
im2.Adv.Fnode = Fnode;
im2.Adv.Fedges = Fedges;
im2.Adv.nodeIn = nodeIn;
im2.Adv.volg = volg;

% store vessel wall flux info
im2.VesFlux.Lflux = Lflux;
im2.VesFlux.gfMap = gfMap;

%should save at this point
save([savefolder,'mesh.mat'],'im2');

end

