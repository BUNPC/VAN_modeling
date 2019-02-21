% This is a revised version of the original code from David.
% 3/12/2012 by L. Gagnon
%
% In this new version I save each frame as a tiff instead of creating a
% movie. It is then easy to reconstruct the movie at the end 6/16/2012.
%
% (from David)
% TO DO
% Implement velocity BC. If velocity not specified then use Lipowsky's data
% on velocity versus diameter for Arterioles and venules. My polynomial fit
% gives 
%    vel (mm/s) = -5.28e-6 d^4 + 2.11e-4 d^3 + 0.0113 d^2 - 0.405 c + 4.70
% where  d=8 um for capillaries, d>8 is venules
% d<8 is arterioles but the true diameter is d_art = (8-d) + 8
%
% Mr and M are sparse. COuld I define them as sparse initially rather than 
% converting them from huge matrices to sparse?
%
% Allow viscosity to vary in segments by varying Hct in segments. Paper I
% reviewed recently for JCBFM suggest that this does not impact flow
% profile
%
% get rid of options
%
% I probably need to get this working on edges as well as segments since
% for O2 adv & diffusion with flow dynamics I will need flow in and out for
% individual edges
%
% Need to check for flow reversal at boundary nodes since this will mess up
% my O2 adv and diffusion.

function [Ftime,Ptime,segVolTime,Fedges,segNodeMap]=flowCircuitEq_dyn( im, tDyn, seg_dvdt_t, segDilationDynType, lstBvis, climvalueVol, climvalueFlow, climvalueDilation, climvaluePressure, lstSegActive, lstSegPassive, print_flag )

if ~exist('options')
    options = [0 0 0];  % options(1) - map pressure from seg node to graph node
                        % options(2) - map edge flow to node flow
end

if ~exist('lstBvis')
    lstBvis = [];
end
if isempty(lstBvis)
    lstBvis = find(im.Mesh.boundary(:,end)==0);
end

if ~isfield(im,'Rscale')
    Rscale = 1;
else
    Rscale=im.Rscale;
    display(sprintf('An Rscale of %d has been used in the flow computation',Rscale))
end

nodeEdges = im.nodeEdges;
nodePos = im.nodePos;
nodeDiam = im.nodeDiam;
nodeBC = im.nodeBC;
nodeBCType = im.nodeBCType;
nodeType = im.nodeType;

nNodes = size(nodePos,1);
nEdges = size(nodeEdges,1);

nSegs = size(im.segEndNodes,1);
segNodes = setdiff(unique(im.segEndNodes),0);
nSegNodes = length(segNodes);
%nSegNodes = max(im.segEndNodes(:));
[foo,segEndNodes] = ismember(im.segEndNodes,segNodes); %convert segEndNodes to reduced segment space
segLen = im.segLen';
segDiam = max(im.segDiam',5);
segVolo = 3.14159*segLen.*(segDiam/2).^2;
segVol = segVolo;
segBeta = 1*ones(nSegs,1); %must be >=1 to ensure stability, higher means stiffer
segTau_v = 1 * ones(nSegs,1);
%segNodeMap = im.segNodeMap; 
%segNodeMap = 1:nNodes;
segNodeMap = segNodes;
segVesType = im.segVesType;
nodePressure = im.nodePressure;
dP_flag=0;

[nB,im] = nBupdate( im );


% Viscosity is ~ 2cP
% convert from cP to mmHg s
% 1 mmHg = 133.3 Pa
% 1 Pa s = 10 Poise = 1000 cP
% 1 Poise = 100 cP
% THEREFORE 1 cP = 0.01 Poise = 1e-3 Pa s = 7.5e-6 mmHg s


% Flow = Mr * Pressure
% Mr is presently defined as nSegs x nSegNodes (specifically seg end nodes)
% This only gives us flow in a segment. It won't give us in-flow and
% out-flow. For this I need pressure defined at segment half nodes.
% My Mr will be constructed as 2*nSegs x (nSegNodes+nSegs)
%    where we have F = [Fin; Fout] and P = [SegEndNodes; SegMidNodes]




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATING RESISTANCE MATRIX
%

% USE SEGMENTS

%hwait = waitbar(0,'Calculating flow : creating resistance matrix');    
nSSN = nSegs+nSegNodes;
nSegs2 = 2*nSegs;
for iS = 1:nSegs
    p1 = segEndNodes(iS,1);
    p2 = segEndNodes(iS,2);
    lstMr(iS,1) = iS + (p1-1)*nSegs2;
    lstMr(iS,2) = iS + (nSegNodes+iS-1)*nSegs2;
    lstMr(iS,3) = nSegs+iS + (nSegNodes+iS-1)*nSegs2;
    lstMr(iS,4) = nSegs+iS + (p2-1)*nSegs2;
end
%close(hwait)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATE PRESSURE
% implement boundary conditions
% to construct M p = y
% where p is the pressure at each node point
% we implement flow conservation at all node points
% except at end points where we can specify pressure
% or specify velocity
% RECALL
%    F = Mr P
%    Mr will be constructed as 2*nSegs x (nSegNodes+nSegs)
%    where we have F = [Fin; Fout] and P = [SegEndNodes; SegMidNodes]
% SO
%    M P = y
%    M is (nSegNodes+nSegs)x(nSegNodes+nSegs) where the bottom nSegs
%      simply comes from Mr where y is dV/dt
%    y is dV/dt for segments mid nodes when nSegB==2
%         0 for segment end nodes when nSegB>2
%         and BC when nSegB==1

%hwait = waitbar( 0, 'Calculating flow : creating pressure matrix' );
nSegB = zeros(nSegNodes,1);
for ii=1:nSegNodes
    nSegB(ii) = length(find(segEndNodes==ii));
end
maxNsegB = max(nSegB);

% USE SEGMENTS
lstM1a = [];
%lstM1aa = [];
lstM1b = [];
lstM1seg = [];
lstM2 = [];
lstM2_fi = [];
lstM2P1 = []; %for dP BC
lstM2P2 = []; %for dP BC
lstM2Pm = []; %for dP BC
lstM3 = [];
lstM3_fi = [];
lstM4a = [];
lstM4b = [];
lstM4seg = [];
lsty2 = [];
lsty2_fi = [];
lsty2P1 = []; %for dP BC
lsty2P2 = []; %for dP BC
lsty2Pm = []; %for dP BC
lsty3 = [];
lsty3seg = [];
lstyRow2 = [];
lstyRow2_fi = [];
lstyRow3 = [];
lsty4 = [];
lsty4sign = [];
lstyRow4 = [];


% LOOP OVER SEGMENT END NODES
for ii=1:nSegNodes %loop arround end segment nodes
    [lstR, lstC] = find(segEndNodes==ii);
    
    if nSegB(ii)>1   % flow conservation
                     % I leave this as 1 rather than 2 in case some segment
                     % end nodes exist with no bifurcation. dV/dt is
                     % handled in the next loop below over nSegs
        for jj=1:length(lstR)
            p1 = segEndNodes(lstR(jj),lstC(jj)); % this is just ii
%            p2 = segEndNodes(lstR(jj),mod(lstC(jj),2)+1);
            lstM1a(end+1) = ii + (p1-1)*nSSN; %this is p1
%            lstM1aa(end+1) = ii + (jj-1)*nSSN + (p1-1)*nSSN*maxNsegB;
            lstM1b(end+1) = ii + (nSegNodes+lstR(jj)-1)*nSSN; %this is the midpoint of the corresponding segment
            lstM1seg(end+1) = lstR(jj); %this is the segment # (required to update resistance)
        end
        
    elseif nSegB(ii)==1  % use a BC
        p1 = segEndNodes(lstR,lstC); %the End node (this is ii)
        p2 = segEndNodes(lstR,mod(lstC,2)+1); %the other node of the segment
        
        
        if (nodeBCType(segNodeMap(ii))==5) %only applies to specific end points selected by the user
        
            dP_flag=1;
            lstM2P1(end+1) = ii + (ii-1)*nSSN;
            lstM2Pm(end+1) = ii + (nSegNodes+lstR-1)*nSSN; %position of the midpoint, this entry of M will be -1 for P1-Pmid=dP
            lstM2P2(end+1) = ii + (p2-1)*nSSN; %position of the other node, this entry of M will be -1 for P1-P2=dP
            lsty2P1(end+1) = ii; %row of the End node
            lsty2P2(end+1) = p2; %position of the other node
            lsty2Pm(end+1) = nSegNodes+lstR; %position of the midnode
            
            
            %NEED TO USE Pressure BC for the first itteration
            % We will assume it is defined by the user since we set
            % everything to user at the end of imView3d_flowCircuitEq
            lstM2_fi(end+1) = ii + (ii-1)*nSSN;
            lstyRow2_fi(end+1) = ii;
            lsty2_fi(end+1) = segNodeMap(ii);
                %                 M(ii,ii) = 1;
                %                 y(ii) = nodeBC(segNodeMap(ii));
            
            
            
            
         
        
              
        elseif nodeBCType(segNodeMap(ii))==1 % specified by user
            lstM2(end+1) = ii + (ii-1)*nSSN;
            lstyRow2(end+1) = ii;
            lsty2(end+1) = segNodeMap(ii);
            %                 M(ii,ii) = 1;
            %                 y(ii) = nodeBC(segNodeMap(ii));
            
            
        elseif nodeBCType(segNodeMap(ii))==3  % specified from literature based on diameter
            lstM3(end+1) = ii + (ii-1)*nSSN;
            lstyRow3(end+1) = ii;
            lsty3(end+1) = segNodeMap(ii);
            lsty3seg(end+1) = lstR;
            %                M(ii,ii) = 1;
            %                y(ii) = getPressure( segDiam(lstR), segVesType(lstR) );
        

          
        elseif nodeBCType(segNodeMap(ii))==2 %velocity (user defined)               
            lstM4a(end+1) = ii + (ii-1)*nSSN;    
            lstM4b(end+1) = ii + (p2-1)*nSSN;
            lstM4seg(end+1) = lstR; %this is the segment # (required to update resistance)
            lstyRow4(end+1) = ii;
            lsty4(end+1) = segNodeMap(ii); %row of the End node
            if lstC==1
                lsty4sign(end+1) = -1;
            else
                lsty4sign(end+1) = 1;
            end
            %         M(ii,p1) = M(ii,p1) + 1/eR(lstR);
            %         M(ii,p2) = M(ii,p2) - 1/eR(lstR);                
            %         vel = nodeBC(segNodeMap(ii));
            %         flow = vel * 3.14159 * (segDiam(lstR)/2)^2;
            %         if lstC==1
            %           y(ii) = -flow;
            %         else
            %         y(ii) = flow;
            %    end
            
        
     
        else % no BC specified so assume vel = 0
            segNodeMap(ii)
            nodeBCType(segNodeMap(ii))
            error( 'we should not get here if only pressure BC employed' )
%             M(ii,p1) = M(ii,p1) + 1/eR(lstR);
%             M(ii,p2) = M(ii,p2) - 1/eR(lstR);            
%             vel=0;            
%             if lstC==1
%                 y(ii) = -vel;
%             else
%                 y(ii) = vel;
%             end
        end        
    end  % End of use a BC
end % End of loop over nodes


yRow3=[]; %for cases where BCs are user defined
for ii=1:length(lstyRow3)  
    yRow3(ii) = getPressure(segDiam(lsty3seg(ii)),floor(segVesType(lsty3seg(ii))));
end
yRow4=[]; %for cases where velocities are specified
for ii=1:length(lstyRow4)   % Maybe I want this fixed during iterations
    yRow4(ii) = lsty4sign(ii) * 1000 * nodeBC(lsty4(ii)) * 3.14159 * (segDiam(lstM4seg(ii))/2)^2; %1000 is to convert mm/s to um/s    
end


% % get list of end segs
% lst = find(nSegB==1);
% lstEndSegA = [];
% lstEndSegC = [];
% lstEndSegV = [];
% for ii=1:length(lst)
%     [lstR,lstC] = find(segEndNodes==lst(ii));
%     if floor(segVesType(lstR))==1
%         lstEndSegA(end+1) = lstR;
%     elseif segVesType(lstR)==2
%         lstEndSegC(end+1) = lstR;
%     elseif segVesType(lst)==3
%         lstEndSegV(end+1) = lstR;
%     end
% end

% find end seg Art with largest diam
%[foo,iAmax]=max(segDiam(lstEndSegA))
%iAmax = lstEndSegA(iAmax);
%lstDilateArt = find(segVesType==1.1);
%lstDilateArt = lstDilateArt(3); %3,4
% 
% lstSegActive = find(segDilationDynType~=1);
% lstSegPassive = find(segDilationDynType==1);



% SOLVE

% initialize for first step
% R = 8 viscosity * len / (pi radius^4)
% divide segLen by 2 to split resistance between in-flow and out-flow

Rsego = Rscale*(128 * 2*7.5e-6 * (segLen/2) ./ (3.14159 * segDiam.^4));
Rseg = Rsego;
segDiamo = segDiam;

Mr = sparse(nSegs*2, nSegNodes+nSegs);
M = sparse(nSSN,nSSN);
y = zeros(nSSN,1);
Po = zeros(nSSN,1);

seg_dVdt = zeros(nSegs,1);

% iterate
dt = tDyn(2)-tDyn(1);
nT = length(tDyn);
Ptime = zeros(nSSN,nT);
VelTime = zeros(nSegs,nT);
Ftime = zeros(nSegs*2,nT);
segVolTime = zeros(nSegs,nT);

%edges
Fedges = zeros(size(im.nodeEdges,1),nT);

%construct arterial dilation vector (for display only)
dilTrace=zeros(size(seg_dvdt_t(:,2:end)));
for iT = 2:nT
    dilTrace(iT,:)=dilTrace(iT-1,:)+seg_dvdt_t(iT,2:end).*dt;
end
dilTrace=sqrt(1+dilTrace)-1; %because we use dV/dt instead of dDiam/dt



for iT = 1:nT
    
    
    disp( sprintf('Iter %d of %d  (%.3f sec)',iT,nT,tDyn(iT)) )
    drawnow
    Mr(lstMr(:)) = [1./Rseg; -1./Rseg; 1./Rseg; -1./Rseg];
    M = sparse(nSSN,nSSN);
    for ii=1:length(lstM1a) % This can be vectorized
        M(lstM1a(ii)) = M(lstM1a(ii)) + 1/Rseg(lstM1seg(ii)); %branching point
        M(lstM1b(ii)) = M(lstM1b(ii)) - 1/Rseg(lstM1seg(ii)); %midpoint
    end
    for ii=1:length(lstM4a) % This can be vectorized
        M(lstM4a(ii)) = M(lstM4a(ii)) + 1/Rseg(lstM4seg(ii)); %branching point
        M(lstM4b(ii)) = M(lstM4b(ii)) - 1/Rseg(lstM4seg(ii)); %midpoint
    end
    
    if (dP_flag)
        if (iT==1)
            M(lstM2_fi) = 1;
            M(lstM3_fi) = 1;
            y(lstyRow2_fi) = nodeBC(lsty2_fi); %these stay fix during itteration
        else
            M(lstM2P1) = 1;
            M(lstM2P2) = -1;
            y(lsty2P1) = DeltaP; 
        end
    end
    
    
    M(lstM2) = 1;
    M(lstM3) = 1;
    y(lstyRow2) = nodeBC(lsty2); %these stay fix during itteration
    y(lstyRow3) = yRow3; %these stay fix during itteration
    y(lstyRow4) = yRow4; %these stay fix during itteration
    
    
    M(nSegNodes+[1:nSegs],:) = Mr(1:nSegs,:)-Mr(nSegs+[1:nSegs],:);
    y(nSegNodes+[1:nSegs]) = seg_dVdt;    
    
    
    scl=sum(M.^2,2).^0.5;
    M = spdiags(1./scl,0,length(scl),length(scl)) * M;
    y = y ./ scl;
%    tic 
%    if iT<2
        P = (M'*M)\M'*y;
%    else
%        [P,flag,relres,iter,resvec] = qmr(M'*M,M'*y,1e-12,1000,[],[],Po);
%        [P,flag,relres,iter,resvec] = bicg(M'*M,M'*y,1e-12,1000,[],[],Po);
%        [P,flag,relres,iter,resvec] = bicgstab(M'*M,M'*y,1e-12,1000,[],[],Po);
%        [P,flag,relres,iter,resvec] = gmres(M'*M,M'*y,[],1e-12,1000,[],[],Po);
%    end
%    toc

    if iT==1
        segAo = (segVolo.^segBeta) ./ P(nSegNodes+[1:nSegs]);
        diffTrace=0;
        DeltaP = P(lsty2P1)-P(lsty2P2);
        nodePressureFlag=0; %this is just for the display
    end
    
    % update arteriole dVdt separately for passive segments
    
    % passive segments
    foo = ( (segAo.*P(nSegNodes+[1:nSegs])).^(1./segBeta) - ...
                 segVol )./segTau_v;
    seg_dVdt(lstSegPassive) = foo(lstSegPassive);
    
    % active segments   
    foo = seg_dvdt_t(iT,:)';
    seg_dVdt(lstSegActive) = segVolo(lstSegActive).*foo(segDilationDynType(lstSegActive));

    
    %     seg_dVdt(lstSegA) = 0;
%     if iT<20
%         seg_dVdt(lstDilateArt) = 1e-1 * segVol(lstDilateArt);
%     else
%         seg_dVdt(lstDilateArt) = 0;
%     end
    
    segVol = segVol + seg_dVdt*dt;
    Rseg = Rsego .* (segVolo./segVol).^2;
    segDiam = segDiamo .* (segVol./segVolo).^(0.5);
    Po = P;
    Ptime(:,iT) = P;
    segVolTime(:,iT) = segVol';
    Ftime(:,iT) = Mr * P;
    
    
    %%%%% Map segment to edges %%%%%%%%%%%%
    % Interpolate Fedges from Flow_In to Flow_Out along segments

    Fseg_In = Ftime(1:nSegs,iT);
    Fseg_Out = Ftime(nSegs+[1:nSegs],iT);
    for iSeg = 1:nSegs
        lstE = find(im.edgeSegN==iSeg);
        NlstE = length(lstE);
        n1 = im.segEndNodes(iSeg,1);
        n2 = im.segEndNodes(iSeg,2);
        nn = n1;
        nold = [];
        eold = [];
        cpt=0;
        while nn~=n2
            cpt=cpt+1;
            lstE1 = find(im.nodeEdges(lstE,1)==nn | im.nodeEdges(lstE,2)==nn);
            eIdx = setdiff( lstE(lstE1), eold );
            eold(end+1) = eIdx;
            if NlstE>1
                Finterp = (cpt-1)/(NlstE-1)*(Fseg_Out(iSeg)-Fseg_In(iSeg))+Fseg_In(iSeg);
            else
                Finterp = (Fseg_Out(iSeg)+Fseg_In(iSeg))./2; %if only one edge then mean flow (not perfect)
            end
            if im.nodeEdges(eIdx,1)==nn
                Fedges(eIdx,iT) = Finterp;
            else
                Fedges(eIdx,iT) = -Finterp;
            end
            
            lstN = im.nodeEdges(lstE(lstE1),:);
            nold(end+1) = nn;
            nn = setdiff( lstN(:), nold );
        end
    end
    

    
    
    
    
    %%%%% FROM HERE IT IS JUST DISPLAY %%%%%%%%%%%%
    
    %flag to display if we want
    if (print_flag && mod(iT,print_flag)==0)
        
        if nodePressureFlag==0
            nodePressureFlag=1;
            %interpolate Pressure to node
            nodePressure1 = zeros(nNodes,1);
            for iSeg = 1:nSegs
                lstE = find(im.edgeSegN==iSeg);
                i1 = find(segNodeMap==im.segEndNodes(iSeg,1));
                i2 = find(segNodeMap==im.segEndNodes(iSeg,2));
                dPint = Ptime(i2,1) - Ptime(i1,1);
                nSteps = length(lstE);
                dPstep = dPint/nSteps;
                n1 = segNodeMap(i1);
                n2 = segNodeMap(i2);
                nodePressure1(n1) = Ptime(i1,1);
                nn = n1;
                Po = P(i1,1);
                nold = [];
                while nn~=n2
                    lstE1 = find(im.nodeEdges(lstE,1)==nn | im.nodeEdges(lstE,2)==nn);
                    lstN = im.nodeEdges(lstE(lstE1),:);
                    nold(end+1) = nn;
                    nn = setdiff( lstN(:), nold );
                    Po = Po + dPstep;
                    nodePressure1(nn) = Po;
                end
            end
        end
        
        %interpolate Pressure to node
        nodePressure = zeros(nNodes,1);
        for iSeg = 1:nSegs
            lstE = find(im.edgeSegN==iSeg);
            i1 = find(segNodeMap==im.segEndNodes(iSeg,1));
            i2 = find(segNodeMap==im.segEndNodes(iSeg,2));
            dPint = P(i2) - P(i1);
            nSteps = length(lstE);
            dPstep = dPint/nSteps;
            n1 = segNodeMap(i1);
            n2 = segNodeMap(i2);
            nodePressure(n1) = P(i1);
            nn = n1;
            Po = P(i1);
            nold = [];
            while nn~=n2
                lstE1 = find(im.nodeEdges(lstE,1)==nn | im.nodeEdges(lstE,2)==nn);
                lstN = im.nodeEdges(lstE(lstE1),:);
                nold(end+1) = nn;
                nn = setdiff( lstN(:), nold );
                Po = Po + dPstep;
                nodePressure(nn) = Po;
            end
        end
        
        %Pressure
%         Pnode2(find(im.nodeSegN~=0))=Ptime(nSegNodes+im.nodeSegN(find(im.nodeSegN~=0)),iT);
%         Pnode1(find(im.nodeSegN~=0))=Ptime(nSegNodes+im.nodeSegN(find(im.nodeSegN~=0)),1);
%         DeltaPnode(find(im.nodeSegN~=0))=Pnode2(find(im.nodeSegN~=0))-Pnode1(find(im.nodeSegN~=0));    
%         PressuretoDisp=Pnode2;
        PressuretoDisp = min(max(nodePressure,climvaluePressure(1)),climvaluePressure(2)); 
        %PressuretoDisp=nodePressure-nodePressure1; 
        
        %Flow
        Fseg2 = ( Ftime(1:nSegs,iT)+Ftime((nSegs+1):(2*nSegs),iT) )./2; %mean of in and out flow
        Fnode2(find(im.nodeSegN~=0)) = abs(Fseg2(im.nodeSegN(find(im.nodeSegN~=0)))); 
        Fseg1 = ( Ftime(1:nSegs,1)+Ftime((nSegs+1):(2*nSegs),1) )./2; %mean of in and out flow;
        Fnode1(find(im.nodeSegN~=0)) = abs(Fseg1(im.nodeSegN(find(im.nodeSegN~=0))));    
        FlowtoDisp = min(max( 100*(abs(Fnode2) - abs(Fnode1))./(abs(Fnode1)),climvalueFlow(1)),climvalueFlow(2));
        %FlowtoDisp = sign(Fnode2-Fnode1).*log10(abs(Fnode2-Fnode1));
        
        %%%%% Uncomment this to plot absolute velocity (in mm/s) instead of
        %%%%% relative flow changes
        %Velocity
%         Velseg=Fseg2./( pi.*(segDiam./2).^2 )/1000; %divide by 1000 to convert um/s to mm/s
%         Velnode(find(im.nodeSegN~=0)) = abs(Velseg(im.nodeSegN(find(im.nodeSegN~=0)))); 
%         FlowtoDisp = min(max(abs(Velnode),climvalueFlow(1)),climvalueFlow(2));
%         
       
        %Volume
        DeltaVol(find(im.nodeSegN~=0))=100*(segVolTime(im.nodeSegN(find(im.nodeSegN~=0)),iT) - ...
            segVolTime(im.nodeSegN(find(im.nodeSegN~=0)),1)) ./ ...
            segVolTime(im.nodeSegN(find(im.nodeSegN~=0)),1);
        VoltoDisp = min(max(DeltaVol,climvalueVol(1)),climvalueVol(2));
        
        
        
            
%             VelSeg = Fseg2([1:nSegs]) ./ (3.14159*(segDiam/2).^2); % Need to use updated Diameter
%            Fnode(find(im.nodeSegN~=0)) = Fseg(im.nodeSegN(find(im.nodeSegN~=0)))
%             VelNode2(find(im.nodeSegN~=0)) = VelSeg(im.nodeSegN(find(im.nodeSegN~=0)));
%              Pnode = VelNode;
%            max(VelSeg)
%             VelTime(:,iT) = VelSeg;

%             Fseg1 = Ftime(:,1);
%             VelSeg = Fseg1([1:nSegs]) ./ (3.14159*(segDiam/2).^2);
%             VelNode1(find(im.nodeSegN~=0)) = VelSeg(im.nodeSegN(find(im.nodeSegN~=0)));
%             Fnode1(find(im.nodeSegN~=0)) = abs(Fseg1(im.nodeSegN(find(im.nodeSegN~=0))));
%             
%            Pnode = min(abs(VelNode1)/1e3,10);
%            Pnode = min((abs(VelNode2) - abs(VelNode1))/1e3,5);
%            Pnode = min(max( (abs(VelNode2) - abs(VelNode1))./abs(VelNode1),-0.5),0.5)
%            Pnode = min((Fnode2 - Fnode1)/1e6,1);
%            Pnode = min((abs(Fnode2)/1e6),15);

%         
            
        
        
        fooFlow = zeros(size(im.Mesh.node,1),1);
        fooVol = zeros(size(im.Mesh.node,1),1);
        fooPressure = zeros(size(im.Mesh.node,1),1);
        lst = find(im.VesFlux.gfMap>0);
        fooFlow(im.Mesh.vesWallnode) = FlowtoDisp(im.VesFlux.gfMap(lst));
        fooVol(im.Mesh.vesWallnode) = VoltoDisp(im.VesFlux.gfMap(lst));
        fooPressure(im.Mesh.vesWallnode) = PressuretoDisp(im.VesFlux.gfMap(lst));
        
       
        if iT/print_flag==1 
            
            hf=figure;
            clf
            set(gcf,'renderer','zbuffer')
            dim_fig=get(hf,'Position');
            dim_fig(3)=1.5*dim_fig(3);
            dim_fig(4)=2.0*dim_fig(4);
            set(hf,'Position',dim_fig);
            %TimeCourseMovie=moviein(nT,hf);
        end
        
        subplot(2,2,1)
            plot(tDyn,100*dilTrace,'b','linewidth',2);hold on; %diam ~ Vol^(1/2)
            yr=climvalueDilation;
            plot([tDyn(iT) tDyn(iT)],[yr(1) yr(2)],'r','linewidth',2)
            ylim(yr);
            xlim([tDyn(1) tDyn(end)])
            xlabel('time (s)','fontsize',14,'fontweight','bold')
            ylabel('arterial dilation (%)','fontsize',14,'fontweight','bold')
            title('Arterial dilation (%)','fontsize',14,'fontweight','bold')
            hold off;
            
        subplot(2,2,2)
            hsurf = trisurf( im.Mesh.boundary(lstBvis,1:3), im.Mesh.node(:,1), im.Mesh.node(:,2), -im.Mesh.node(:,3), fooFlow, 'linestyle','none','facecolor','flat');
            set(gca,'clim',climvalueFlow)
            colorbar
            axis image
            grid off
            colormap(jet(64))
            view(0,90)            
            title('Flow change (%)','fontsize',14,'fontweight','bold')
%             title('Velocity (mm/s)','fontsize',14,'fontweight','bold')
        subplot(2,2,3)
            hsurf = trisurf( im.Mesh.boundary(lstBvis,1:3), im.Mesh.node(:,1), im.Mesh.node(:,2), -im.Mesh.node(:,3), fooPressure, 'linestyle','none','facecolor','flat');
            set(gca,'clim',climvaluePressure)
            colorbar
            axis image
            grid off
            colormap(jet(64))
            view(0,90)            
            title('Pressure (mmHg)','fontsize',14,'fontweight','bold')
        subplot(2,2,4)
            hsurf = trisurf( im.Mesh.boundary(lstBvis,1:3), im.Mesh.node(:,1), im.Mesh.node(:,2), -im.Mesh.node(:,3), fooVol, 'linestyle','none','facecolor','flat');
            set(gca,'clim',climvalueVol)
            colorbar
            axis image
            grid off
            colormap(jet(64))
            view(0,90)            
            title('Volume change (%)','fontsize',14,'fontweight','bold')
        
        %print the frame to create a movie at the end
        set(gcf,'PaperPositionMode','auto')
        eval(sprintf('print -dtiff -r0 frames_flow/FlowVolume_%d.tiff',iT))
        
    end
        
        if iT>1
            pause(0.1)
        end
end

return





% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Calculate Velocity
% Fseg = Mr * P;
% 
% % map flow in segments to flow in edges
% % make sure we get the sign right
% Fsign = sign(im.nodePressure(im.nodeEdges(:,1)) - im.nodePressure(im.nodeEdges(:,2)))';
% F = Fsign .* abs(Fseg(im.edgeSegN));
% Fnode(find(im.nodeSegN~=0)) = Fseg(im.nodeSegN(find(im.nodeSegN~=0)))
% nodeVel = Fnode ./ (3.14159*(im.segDiam/2).^2);
% 
% 
% 
% im.segFlow = Fseg;
% im.segVel = Fseg ./ (3.14159*(im.segDiam(:)/2).^2);
% im.edgeFlow = F;
% im.edgeVel = im.segVel(im.edgeSegN);
% 
% 
% 
% % This is wrong, check imView3d_flowCircuitEq
% if im.flagUseSegments
%     Fseg = F;
%     F = Fseg(im.edgeSegN);
%     
%     im.segFlow = Fseg;
%     im.segVel = Fseg ./ (3.14159*(im.segDiam(:)/2).^2);
%     im.edgeFlow = F;
%     im.edgeVel = im.segVel(im.edgeSegN);
% end
% Fedges = F;
% 
% 
% return
% 
% 
% 
% 
% Fnode = zeros(nNodes,1);
% if options(2)
%     hwait = waitbar(0,'Calculating flow : mapping edge flow to node flow');
%     for ii=1:nNodes
%         waitbar(ii/nNodes,hwait);
%         [lstR,lstC] = find(nodeEdges==ii);
%         %    Fnode(ii) = mean(F(lstC).*((-1).^(lstC+1)));
%         Fnode(ii) = mean(abs(F(lstR)));
%         %    F(lstR)'
%         %    nodeEdges(lstR,:)
%         %    pause
%     end
%     close(hwait)
% end
% 
% 
% if ~im.flagUseSegments
% 
%     % USE EDGES
% 
%     im.nodeVel = zeros(nNodes,1);
%     hwait = waitbar(0,'Calculating flow : calculate velocity step 2 of 3');
%     for iN=1:nNodes
%         waitbar(iN/nNodes,hwait);
%         im.nodeVel(iN) = Fnode(iN) / (3.14159*nodeDiam(iN)^2/4);
%     end
%     close(hwait)
%     im.edgeFlow = Fedges;
%     im.edgeVel = zeros(nEdges,1);
%     hwait = waitbar(0,'Calculating flow : calculate velocity step 3 of 3');
%     for iE=1:size(nodeEdges,1)
%         waitbar(iE/nEdges,hwait);
%         rad = mean(nodeDiam(nodeEdges(iE,:)))/2;
%         im.edgeVel(iE) = Fedges(iE) / (3.14159*rad^2);
%     end
%     close(hwait)
% 
% else
% 
%     % USE SEGMENTS
%     
%     im.nodeVel = zeros(nNodes,1);
%     hwait = waitbar(0,'Calculating flow : calculate velocity step 2 of 3');
%     for iN=1:nNodes
%         waitbar(iN/nNodes,hwait);
%         if nB(iN)==2
%             im.nodeVel(iN) = Fnode(iN) / (3.14159*segDiam(im.nodeSegN(iN))^2/4);
%         else
%             im.nodeVel(iN) = 0;
%         end
%     end
%     close(hwait)
%     im.edgeFlow = Fedges;
%     im.edgeVel = zeros(nEdges,1);
%     hwait = waitbar(0,'Calculating flow : calculate velocity step 3 of 3');
%     for iE=1:size(nodeEdges,1)
%         waitbar(iE/nEdges,hwait);
% %        rad = mean(segDiam(im.nodeSegN(nodeEdges(iE,:))))/2;
%         rad = max(segDiam(im.nodeSegN(nodeEdges(iE,:))))/2;
%         im.edgeVel(iE) = Fedges(iE) / (3.14159*rad^2);
%     end
%     close(hwait)
% 
% end
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % quiver plot of flow
% %%
% % xq = [];
% % yq = [];
% % zq = [];
% % uq = [];
% % vq = [];
% % wq = [];
% % for ii=1:size(nodeEdges,1)
% %     xq(end+1) = mean(nodePos(nodeEdges(ii,:),1));
% %     yq(end+1) = mean(nodePos(nodeEdges(ii,:),2));
% %     zq(end+1) = mean(nodePos(nodeEdges(ii,:),3));
% %     rq = nodePos(nodeEdges(ii,2),:) - nodePos(nodeEdges(ii,1),:);
% %     rq = rq / norm(rq);
% %     uq(end+1) = F(ii)*rq(1);
% %     vq(end+1) = F(ii)*rq(2);
% %     wq(end+1) = F(ii)*rq(3);
% % end
% % figure(10);
% % quiver3(xq,yq,zq,uq,vq,wq,1,'linewidth',1.5,'maxheadsize',1)
% % 
% 
% 
