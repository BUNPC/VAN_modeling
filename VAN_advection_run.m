function [c,cg,cbg,t,AmatFlux,BmatFlux, OCoutput] = VAN_advection_run( im, nIter, nIterPerPlot, AmatFlux, BmatFlux, c, cg, cbg, tStart, BCvec )

% To Do
%
% update Ctg, put the 40 as a param in im.
%
% what should vessel wall permeability be?
%
% setup the simVessel script to save the necesary mesh info so that it is
%    only generated once. This will allow me to easily change velocity and
%    not regenerate the mesh.... but Ladv needs to be regenerated.  Hmmm...
%    so I need a function that generates Ladv based on changes in flow and
%    volume along the graph. I will also need to change Lflux to account
%    for diameter changes that will modify the surface area of the flux.
%    Fortunately, I won't need to change the oxygen diffusion matrices
%
% so2_func is called a lot. Can we make this faster? Use profile to see if
%    it is expensive
%
% Diffusion - Modify the permeability coefficient to handle vessel dilation
%   without having to redo the mesh at each itteration (8/16/2012 LG).
%
% Advection - put a check in for the stability criteria
%


if ~isempty(cbg)
    flagInitC = 0;
    if size(cbg,2)>1 | size(cg,2)>1 | size(c,2)>1 | size(cbg,1)~=size(cg,1) 
        error('USAGE: If passing c, cg, and cbg, then\nthey must be column vectors and length(cg)=length(cbg)')        
    end
else
    flagInitC = 1;
end

%add that to avoid error LG 1/10/2012
if ~exist('tStart')
    tStart=0;
end

if ~isfield(im,'nIterPerTrecord')
    nIterPerTrecord = 0;
else
    nIterPerTrecord = im.nIterPerTrecord;
end
if nIterPerTrecord==0
    nIterPerTrecord=nIter*2;
end

if ~isfield(im,'nIterPerOCcheck')
    nIterPerOCcheck = 1;
else
    nIterPerOCcheck = im.nIterPerOCcheck;
end

nAdvIterPerIter = im.nAdvIterPerIter;

if ~exist('AmatFlux')
    AmatFlux = [];
end
if ~exist('BmatFlux')
    BmatFlux = [];
end

if ~isfield(im,'species')
    species=1;
else
    species=im.species;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract values of interest
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nng = length(im.Adv.volg);
nn = size(im.Mesh.nodevol,1)-nng;

alpha = im.alpha;
% po2a = im.po2a;
po2t = im.po2t;
OCo = im.OCo;
Chb = im.Chb;
Hct = im.Hct;
hwall = im.hwall;
hwallalpha = im.hwall * im.alpha;
Kves = im.Kves;
KvesArtFact = im.KvesArtFact;

dt = im.dt;
Do2 = im.Do2;
dcoeff=Do2*ones(nn,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up FEM matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sita=0.5;

[Amat,Bmat]=diffusion2(dt,sita,im.Mesh.node,im.Mesh.elem,im.Mesh.iiPrepMesh,...
    im.Mesh.jjPrepMesh,im.Mesh.connnum,im.Mesh.vol,dcoeff,nng);

if (isempty(AmatFlux) || isempty(BmatFlux))
    disp( 'Creating Flux Matrices' )
    tic
    % C^{k+1} - C^k = dt (area/Vol) (-Kves( C_out - C_in ) / (h alpha))
    k  = im.Mesh.nodearea(im.Mesh.vesWallnode)*Kves*dt/hwallalpha;
    
    %different permeability for arteries
    type1WallNodes = im.Mesh.vesWallnode;
    segVesType = im.segVesType';
    lstNart = find( segVesType(im.nodeSegN(im.VesFlux.gfMap(type1WallNodes)))'==1 );
    k(lstNart) = KvesArtFact * im.Mesh.nodearea(lstNart)*Kves*dt/hwallalpha;
    
    kt = k;%./(nodevol(vesWallnode));
    kv = k;%./volg(gfMap(vesWallnode));
    v1 = im.Mesh.nodevol(im.Mesh.vesWallnode);
    v2 = im.Adv.volg;%(gfMap(vesWallnode));
    %v1 = ones(size(v1));
    %v2 = ones(size(v2));
    %nPatchesPerVesselNode = sum(Lflux,2);
    kvg = kv;% ./ nPatchesPerVesselNode(gfMap(vesWallnode));
    AmatFlux = sparse(nn+nng,nn+nng);
    BmatFlux = sparse(nn+nng,nn+nng);
    
    lst1 = sub2ind(size(AmatFlux),im.Mesh.vesWallnode,im.Mesh.vesWallnode);
    lst2 = sub2ind(size(AmatFlux),im.Mesh.vesWallnode,nn+im.VesFlux.gfMap(im.Mesh.vesWallnode));
    lst3 = sub2ind(size(AmatFlux),nn+im.VesFlux.gfMap(im.Mesh.vesWallnode),im.Mesh.vesWallnode);
    lst4 = sub2ind(size(AmatFlux),nn+[1:nng],nn+[1:nng]);
    AmatFlux(lst1) = v1 + sita*kt;
    AmatFlux(lst2) = -sita*kt;
    AmatFlux(lst3) = -sita*kv;
    AmatFlux(lst4) = v2+sita*(im.VesFlux.Lflux*kvg);
    BmatFlux(lst1) = v1 - (1-sita)*kt;
    BmatFlux(lst2) = (1-sita)*kt;
    BmatFlux(lst3) = (1-sita)*kv;
    BmatFlux(lst4) = v2-(1-sita)*(im.VesFlux.Lflux*kvg);
    
    
    
    % Input Nodes: type 1 boundary condition
    lstIn = find(im.Adv.nodeIn==1);
    for ii=1:length(lstIn)
        AmatFlux(nn+lstIn(ii),1:nn) = 0;
        AmatFlux(1:nn,nn+lstIn(ii)) = 0;
        AmatFlux(nn+lstIn(ii),nn+lstIn(ii)) = 1;
        lstSurfNodes = find(im.VesFlux.gfMap==lstIn(ii));
        lst = sub2ind(size(AmatFlux),lstSurfNodes,lstSurfNodes);
        AmatFlux(lst) = 1;
        
        BmatFlux(nn+lstIn(ii),1:nn) = 0;
        BmatFlux(1:nn,nn+lstIn(ii)) = 0;
        BmatFlux(nn+lstIn(ii),nn+lstIn(ii)) = 1;
        lstSurfNodes = find(im.VesFlux.gfMap==lstIn(ii));
        lst = sub2ind(size(BmatFlux),lstSurfNodes,lstSurfNodes);
        BmatFlux(lst) = 1;
    end
    % Output Nodes
    lstOut = find(im.Adv.nodeIn==-1);
    for ii=1:length(lstOut)
        AmatFlux(nn+lstOut(ii),1:nn) = 0;
        AmatFlux(1:nn,nn+lstOut(ii)) = 0;
        AmatFlux(nn+lstOut(ii),nn+lstOut(ii)) = 1;
        lstSurfNodes = find(im.VesFlux.gfMap==lstOut(ii));
        lst = sub2ind(size(AmatFlux),lstSurfNodes,lstSurfNodes);
        AmatFlux(lst) = 1;
        
        BmatFlux(nn+lstOut(ii),1:nn) = 0;
        BmatFlux(1:nn,nn+lstOut(ii)) = 0;
        BmatFlux(nn+lstOut(ii),nn+lstOut(ii)) = 1;
        lstSurfNodes = find(im.VesFlux.gfMap==lstOut(ii));
        lst = sub2ind(size(BmatFlux),lstSurfNodes,lstSurfNodes);
        BmatFlux(lst) = 1;
    end
    toc
else
    disp( 'Using passed Flux Matrices' );
end

disp( 'Adding All Matrices together' )
tic
AmatAll = sparse((nn+nng),nn+nng);
BmatAll = sparse((nn+nng),nn+nng);
AmatAll = [Amat + AmatFlux];
BmatAll = [Bmat + BmatFlux];
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SET INITIAL CONDITIONS (modified by L. Gagnon)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp( 'Setting initial conditions' )
tic


VlstBCIn = BCvec(:,1);
po2_in = BCvec(:,2);
VlstBCIn_all = VlstBCIn; %this is for display only
VlstBCOut = find(im.Adv.nodeIn==-1);
VlstNoBC = find(im.Adv.nodeIn~=1);
if flagInitC %initialize matrices and apply BCs
    disp( sprintf('Fixing pO2 input on %d nodes',length(VlstBCIn)) )
    c=ones(nn,1)*po2t * alpha;
    cg = zeros(nng,1);
    cbg = zeros(nng,1);
    
    %set concentration on all vessels based on tissue pO2
    cg(:) = po2t*alpha;
    cbg(:) = 4*Chb*Hct*so2_func(cg/alpha,species);
        
    %force pO2 for vessels for which we have a BC
    cg(VlstBCIn) = po2_in*alpha; 
    cbg(VlstBCIn) = 4*Chb*Hct*so2_func(cg(VlstBCIn)/alpha,species);
    
   



else %use the matrices we passed as initial conditions but we are still forcing BCs for inflow nodes
    disp( 'pO2 input determined from previous c, cg, cbg' )
    %cg(VlstBCIn) = po2a*alpha;
    cg(VlstBCIn) = po2_in*alpha;
    cbg(VlstBCIn) = 4*Chb*Hct*so2_func(cg(VlstBCIn)/alpha,species);
end
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RUN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
disp( 'Running...' )

iOCtmp = find(im.tOC<=(1*dt+tStart),1,'last');
os=zeros(size(AmatAll,1),1);
oc = zeros(nn+nng,1);
if size(OCo,1)==1
    oc(1:nn) = OCo(:,iOCtmp) .* im.Mesh.nodevol(1:nn) .* c ./(c+1e-15); %input homogeneous CMRO2
    disp('Initial CMRO2 was homogeneous')
else
    oc(1:nn) = OCo(:,iOCtmp) .* im.Mesh.nodevol(1:nn); %input heterogeneous CMRO2, no roll-off
    disp('Initial CMRO2 was heterogeneous')    
end
iOC = 1;
iOCcheck = 1;
OCoutput=zeros(nIter,1); %total oxygen consummed at each time step


%interpolate flow and volume changes over entire time
tic
disp('Interpolating Flow and Volume changes')
tFlowInt=(im.tFlow(1):dt:im.tFlow(end))'; %this should give something with length of nIter+1
FedgesInt=interp1(im.tFlow',im.Fedges',tFlowInt')';
segRelVolInt=interp1(im.tFlow',im.segRelVol',tFlowInt')';
FtimeInt=interp1(im.tFlow',im.Ftime',tFlowInt')';


%Make sure we don't get NaN because of different time courses
if ~isempty(find(isnan(segRelVolInt))~=0) || ~isempty(find(isnan(FedgesInt))~=0) 
    error('Interpolation resulted in some NaN quantitites. Most probable cause is nIter exceed Flow time course')    
end

%Make sure we have nIter+1 flow and volume values (remember A(k+1) and B(K))
if size(FedgesInt,2)<(nIter+1)
    error('The interpolated Flow time course is not length nIter+1. Please set nIter to be length(tFlowInt)-1')
end    
toc


cgo = cg;
cbgo = cbg;
co = c;
csurfo = [];
csurf = [];
for ii=1:nng;
    lst = find(im.VesFlux.gfMap(im.Mesh.vesWallnode)==ii);
    csurfo(ii)=sum(co(lst).*im.Mesh.nodevol(im.Mesh.vesWallnode(lst)))/sum(im.Mesh.nodevol(im.Mesh.vesWallnode(lst)));
end

iTrecord = 0;
nTrecord = floor(nIter / nIterPerTrecord);
if nTrecord>0
    cTrec = zeros(length(c),nTrecord);
    cgTrec = zeros(length(cg),nTrecord);
    cbgTrec = zeros(length(cbg),nTrecord);
end
toc

for i=1:nIter
    tic
    fprintf(1,'solving iter %d ...\n',i);
    
    %Hold constant concentration for vessels for which we have a BC,
    %otherwise concentration will decrease just by increasing volume. LG 6/03/2012
    cg(VlstBCIn) = po2_in*alpha; 
    cbg(VlstBCIn) = 4*Chb*Hct*so2_func(cg(VlstBCIn)/alpha,species);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Do the oxygen diffusion in the tissue, 
    % and flux across the vessel wall
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	cnew=diffusion2solve(dt,AmatAll,BmatAll,oc,os,[c;cg],im.Mesh.nodevol);
    c = cnew(1:nn);
    cg = cnew(nn+[1:nng]);
    OCoutput(i) = sum(oc(1:nn)); %total oxygen consummed

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Do the advection (modified by L. Gagnon 6/12/2012)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    iFlow = find(tFlowInt<=( (i-1)*dt+tStart ),1,'last');
    
    neg = size(im.nodeEdges,1);
    nng = size(im.nodePos,1);
    
    Ladv_A = zeros(nng,nng);
    Ladv_B = zeros(nng,nng);
    
    Fedges_A = FedgesInt(:,iFlow+1); %This is interpolated in flowCircuitEq to handle Flow_In ~= Flow_Out during dilation
    Fedges_A(neg+1) = 0; % use neg+1 index for 0 flow
    Fedges_B = FedgesInt(:,iFlow);
    Fedges_B(neg+1) = 0; % use neg+1 index for 0 flow
    
    %NOTE: If there is only one edge then I need to separate explicitly
    %flow_in and flow out for this segment because Fedges is the mean flow
    %in this case.
    nSeg=size(FtimeInt,1)/2;
    Ftime_In_A = FtimeInt(1:nSeg,iFlow+1);
    Ftime_Out_A = FtimeInt(nSeg+[1:nSeg],iFlow+1);
    Ftime_In_B = FtimeInt(1:nSeg,iFlow);
    Ftime_Out_B = FtimeInt(nSeg+[1:nSeg],iFlow);
    
    %   We don't want the flow reverse direction in the middle on an
    %   itteration. Flow reversal is handeled by recomputing the advection
    %   matrices at each itteration, but this direction must be constant
    %   between A and B.
    Fedges_A(sign(Fedges_A).*sign(Fedges_B)<0)=0;
    Ftime_In_A(sign(Ftime_In_A).*sign(Ftime_In_B)<0)=0;
    Ftime_Out_A(sign(Ftime_Out_A).*sign(Ftime_Out_B)<0)=0;
    
    volg_A = im.Adv.volg.*segRelVolInt(im.nodeSegN,iFlow+1); 
    volg_B = im.Adv.volg.*segRelVolInt(im.nodeSegN,iFlow); 
    
    eg = im.nodeEdges;
    
    
    
    for ii=1:nng %loop on nodes
        [ir,ic] = find( eg==ii );
        for jj=1:length(ir) %loop on edges
            n1 = eg(ir(jj),1);
            n2 = eg(ir(jj),2);
            eIdx = ir(jj);
            if (n2==ii && Fedges_B(eIdx)>=0) %inFlow
                
                if im.segNedges(im.edgeSegN(eIdx))==1 
                    Ladv_A(ii,n1) = abs(Ftime_Out_A(im.edgeSegN(eIdx))) / volg_A(n1);%inflowing node, so use segment flow_out
                    Ladv_B(ii,n1) = abs(Ftime_Out_B(im.edgeSegN(eIdx))) / volg_B(n1);
                else
                    Ladv_A(ii,n1) = abs(Fedges_A(eIdx)) / volg_A(n1);
                    Ladv_B(ii,n1) = abs(Fedges_B(eIdx)) / volg_B(n1);
                end
                
            elseif (n1==ii && Fedges_B(eIdx)<0) %inFlow
                
                if im.segNedges(im.edgeSegN(eIdx))==1
                    Ladv_A(ii,n2) = abs(Ftime_Out_A(im.edgeSegN(eIdx))) / volg_A(n2);%inflowing node, so use segment flow_out
                    Ladv_B(ii,n2) = abs(Ftime_Out_B(im.edgeSegN(eIdx))) / volg_B(n2);
                else
                    Ladv_A(ii,n2) = abs(Fedges_A(eIdx)) / volg_A(n2);
                    Ladv_B(ii,n2) = abs(Fedges_B(eIdx)) / volg_B(n2);
                end
                
            elseif (n2==ii && Fedges_B(eIdx)<0) || (n1==ii && Fedges_B(eIdx)>=0) %outFlow
                if im.segNedges(im.edgeSegN(eIdx))==1
                    Ladv_A(ii,ii) = Ladv_A(ii,ii) -1*abs(Ftime_In_A(im.edgeSegN(eIdx))) / volg_A(ii);%outflowing node, so use segment flow_in
                    Ladv_B(ii,ii) = Ladv_B(ii,ii) -1*abs(Ftime_In_B(im.edgeSegN(eIdx))) / volg_B(ii);
                else
                    Ladv_A(ii,ii) = Ladv_A(ii,ii) -1*abs(Fedges_A(eIdx)) / volg_A(ii);
                    Ladv_B(ii,ii) = Ladv_B(ii,ii) -1*abs(Fedges_B(eIdx)) / volg_B(ii);
                end
                
            end
        end
        if length(ir)==1 && ((n2==ii && Fedges_B(eIdx)>=0) || (n1==ii && Fedges_B(eIdx)<0))%this is an End Node with inFlow
            
            if im.segNedges(im.edgeSegN(eIdx))==1
                Ladv_A(ii,ii) =  -1*abs(Ftime_Out_A(im.edgeSegN(eIdx))) / volg_A(ii);%inflowing node, so use segment flow_out
                Ladv_B(ii,ii) =  -1*abs(Ftime_Out_B(im.edgeSegN(eIdx))) / volg_B(ii);
            else
                Ladv_A(ii,ii) =  -1*abs(Fedges_A(eIdx)) / volg_A(ii);
                Ladv_B(ii,ii) =  -1*abs(Fedges_B(eIdx)) / volg_B(ii);
            end
            
        elseif length(ir)==1 && ((n2==ii && Fedges_B(eIdx)<0) || (n1==ii && Fedges_B(eIdx)>=0))%End node with outFlow
            Ladv_A(ii,ii) = 0;
            Ladv_B(ii,ii) = 0;
            
        end
    end
    Ladv_A = sparse(Ladv_A);
    AmatAdv = sparse(speye(nng,nng) - sita*Ladv_A*dt/nAdvIterPerIter);
    
    Ladv_B = sparse(Ladv_B);
    BmatAdv = sparse(speye(nng,nng) + (1-sita)*Ladv_B*dt/nAdvIterPerIter);
    
    %convert concentration to number of particle for advection
    nbg = cbg .* volg_B;
    ng = cg .* volg_B;    
    
    invAmatAdvBmatAdv=AmatAdv\BmatAdv;
    for iAdv = 1:nAdvIterPerIter
        nbg = invAmatAdvBmatAdv*nbg;
        ng = invAmatAdvBmatAdv*ng;    
    end
    
    %extract concentration from number of particle
    cbg = nbg ./ volg_A;
    cg = ng ./ volg_A;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Constraints
    % c >= 0 in tissue
    % OC goes to zero as c goes to 0
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    c = max(c,0);
    cg = max(cg,0);
    %cbg = max(cbg,0);
    
    % Do infrequently. This one line seemed to cost 20% of the time
    iOCtmp = find(im.tOC<=(i*dt+tStart),1,'last');
    if iOC~=iOCtmp | iOCcheck==nIterPerOCcheck
        iOCcheck = 0;
        iOC = iOCtmp;
        if size(OCo,1)==1
            oc(1:nn) = OCo(:,iOC) .* im.Mesh.nodevol(1:nn) .* c ./(c+1e-15);  % taken from Aubert2002
                        % role off occurs at 1 uM of oxygen. our units are
                        % umole / um^3, so multiple by 1e-15 um^3 / L
                        % this is about 1 mmHg
        else
            oc(1:nn) = OCo(:,iOC) .* im.Mesh.nodevol(1:nn);  % no role off if heterogeneous CMRO2     
        end
    end
    iOCcheck = iOCcheck + 1;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Update vessel Ct
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for ii=1:60
        cbgnew = 4*Chb*Hct*so2_func(cg(VlstNoBC)/alpha,species);
        dcbg = (cbgnew - cbg(VlstNoBC))/40;
        cbg(VlstNoBC) = cbg(VlstNoBC) + dcbg;
        cg(VlstNoBC) = cg(VlstNoBC) - dcbg;
    end
    cg = real(cg); % these are added because so2_func can return imag if po2<0
    cbg = real(cbg);
    
    cg = max(cg,0);
    cbg = max(cbg,0);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Update dynamic PO2 BCs (never used that ???)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if 0
        cg(VlstBCIn_dyn(:,1)) = cg(VlstBCIn_dyn(:,2));
        cbg(VlstBCIn_dyn(:,1)) = cbg(VlstBCIn_dyn(:,2));
    end
    
  
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Record temporal variation if requested
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if mod(i,nIterPerTrecord)==0
        iTrecord = iTrecord + 1;
        cTrec(:,iTrecord) = c;
        cgTrec(:,iTrecord) = cg;
        cbgTrec(:,iTrecord) = cbg;
        t(iTrecord) = i*dt + tStart;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Visualize
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(mod(i,nIterPerPlot)==0) 
        
        % Fnode is funny as calculated on line 90 in setupMesh.m
        % but only for branch nodes
        fprintf(1,'     O2 Consumption - AVb=%.3e   AVt=%.3e   Tis=%.3e\n',...
            4*Chb*Hct*( sum(abs(im.Adv.Fnode(VlstBCIn_all)).*so2_func(cg(VlstBCIn_all)/alpha,species)) - sum(abs(im.Adv.Fnode(VlstBCOut)).*so2_func(cg(VlstBCOut)/alpha,species)) ),...
            sum(abs(im.Adv.Fnode(VlstBCIn_all)).*cg(VlstBCIn_all)) - sum(abs(im.Adv.Fnode(VlstBCOut)).*cg(VlstBCOut)) + 4*Chb*Hct*( sum(abs(im.Adv.Fnode(VlstBCIn_all)).*so2_func(cg(VlstBCIn_all)/alpha,species)) - sum(abs(im.Adv.Fnode(VlstBCOut)).*so2_func(cg(VlstBCOut)/alpha,species)) ),...
            sum(oc(1:nn)) );
       
    end
    toc
end

if iTrecord>0
    c = cTrec;
    cg = cgTrec;
    cbg = cbgTrec;
end

