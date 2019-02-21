% This script computes the flow-volume dynamics using Anna's 2PM data as
% dilation inputs (rat data from Tian 2010). We use the stack with no capillary end segments (NCES)
% and using constant pressure BC's for arteries and veins.
%
% 3/12/2013 by L. Gagnon
%
function VAN_calcFlow_sigmoid(savefolder,mousefolder1,VesselDiPercent)


%load data with mesh (mesh is just to display during computation)
load([savefolder,'mesh.mat']);

%load branching order
load([mousefolder1,'brOrder.mat']);
im2.grpStat=grpStat;
im2.edgeBRorderArt=edgeBRorderArt;
im2.edgeBRorderVeins=edgeBRorderVeins;


nSegs = size(im2.segEndNodes,1);
segNodes = reshape(im2.segEndNodes,[nSegs*2 1]);
segPosX = reshape(im2.Hvox(1).*im2.nodePos(segNodes,1),[nSegs 2]);
segPosY = reshape(im2.Hvox(2).*im2.nodePos(segNodes,2),[nSegs 2]);
segPosZ = reshape(im2.Hvox(3).*im2.nodePos(segNodes,3),[nSegs 2]);
segPosR = sum( (im2.Hvox(1).*im2.nodePos(im2.segEndNodes(:,1),:) - im2.Hvox(2).*im2.nodePos(im2.segEndNodes(:,2),:)).^2,2).^0.5;
segPosR=max(segPosR,0.1);
segPosTheta = acos( abs(diff(segPosZ,1,2)) ./ segPosR ) * 180/3.14159;
segVesType = round(im2.segVesType)';


%%%% Let's divide the active dilation compartments %%%%%%

% I already computed braching order from imView3d so I will use that

% Will keep these the same 
% (1) arteries 
% (2) capillaries 
% (3) veins 

% Will re-divide arteries and pre-capillary arterioles

% (4) Surface pial artery
lstPial = find(segVesType==1 & segPosTheta>45 & mean(segPosZ,2)<200);
segVesType(lstPial) = 4;

% I have data from Anna that distinguish trunk dilation at each depth

% (5) Diving trunk 0 to 150um
lstDiving0_150 = find(im2.grpStat.segBranchOrder(:,2)==0 & segPosTheta<45 & mean(segPosZ,2)<150);
segVesType(lstDiving0_150) = 5;

% (6) Diving trunk >150um
lstDiving150more = find(im2.grpStat.segBranchOrder(:,2)==0 & segPosTheta<45 & mean(segPosZ,2)>=150);
segVesType(lstDiving150more) = 6;

% (12) 1st branching arterioles < 150
lst1BrOrd0_150 = find(im2.grpStat.segBranchOrder(:,2)==1 & mean(segPosZ,2)<150);
segVesType(lst1BrOrd0_150) = 7;

% (13) 1st branching arterioles > 150
lst1BrOrd150more = find(im2.grpStat.segBranchOrder(:,2)==1 & mean(segPosZ,2)>=150);
segVesType(lst1BrOrd150more) = 8;

% (14) 2nd-3rd-4th branching arterioles < 150
lst234BrOrd0_150 = find( (im2.grpStat.segBranchOrder(:,2)==2 | im2.grpStat.segBranchOrder(:,2)==4 | im2.grpStat.segBranchOrder(:,2)==5) &  mean(segPosZ,2)<150);
segVesType(lst234BrOrd0_150) = 9;

% (15) 2nd-3rd-4th branching arterioles > 150
lst234BrOrd150more = find( (im2.grpStat.segBranchOrder(:,2)==2 | im2.grpStat.segBranchOrder(:,2)==4 | im2.grpStat.segBranchOrder(:,2)==5) &  mean(segPosZ,2)>=150);
segVesType(lst234BrOrd150more) = 10;

%%%%%% Re-sample and assign dilation traces %%%%%%%%%
dtFlow=0.025;
maxTFlow=6;
tFlow=(dtFlow:dtFlow:maxTFlow)';

% %load rat data set (for surface, 1st branches, 2-4th branches)
% load anna_MPM_timecourses.mat
% 
% %remove NaN values
% gidx=find(isnan(ysurface)==0 & isnan(y0branch1)==0 & isnan(y0branch234)==0 & isnan(y150branch1)==0 & isnan(y150branch234)==0 & isnan(y0trunk)==0 & isnan(y150trunk)==0);
% x0=x0(gidx);
% ysurface=ysurface(gidx);
% y0branch1=y0branch1(gidx);
% y0branch234=y0branch234(gidx);
% y150branch1=y150branch1(gidx);
% y150branch234=y150branch234(gidx);
% y0trunk=y0trunk(gidx);
% y150trunk=y150trunk(gidx);
% 
% %resample to tFlow
% ysurface=interp1(x0,ysurface,tFlow);
% y0branch1=interp1(x0,y0branch1,tFlow);
% y0branch234=interp1(x0,y0branch234,tFlow);
% y150branch1=interp1(x0,y150branch1,tFlow);
% y150branch234=interp1(x0,y150branch234,tFlow);
% y0trunk=interp1(x0,y0trunk,tFlow);
% y150trunk=interp1(x0,y150trunk,tFlow);

%create sigmoid dilation traces
a1=3;
a2=2;
a3=VesselDiPercent; %percent changes in amplitude
dilation_trace = a3./(1+exp(-a1.*(tFlow-a2)));

%assign active dilation
segDilationDynType = segVesType;
%segDilationDyn = [zeros(length(tFlow),3) ysurface y0trunk y150trunk y0branch1 y150branch1 y0branch234 y150branch234]/100;
segDilationDyn = [zeros(length(tFlow),3) dilation_trace dilation_trace dilation_trace dilation_trace dilation_trace dilation_trace dilation_trace]/100;
%segDilationDyn = filtfilt(ones(10,1),10,segDilationDyn); %LPF to make things smooth
seg_dvdt_t = zeros(size(segDilationDyn));
seg_dvdt_t(2:end,:) = diff((1+segDilationDyn).^2,1,1)/(dtFlow); %square because volume ~ diameter^2
lstSegActive = find(segDilationDynType>=4);
lstSegPassive = find(segDilationDynType<=3);


%%%% Setup mesh for display %%%%%%%%

foo = zeros(size(im2.Mesh.node,1),1);
lst = find(im2.VesFlux.gfMap>0);
foo(im2.Mesh.vesWallnode) = segVesType(im2.nodeSegN(im2.VesFlux.gfMap(lst)))';% Vessel type refine for descending, branching, etc

lstNav = find( segVesType(im2.nodeSegN(im2.VesFlux.gfMap(lst)))'~=2 ); %arteries and veins
lstNc = find( segVesType(im2.nodeSegN(im2.VesFlux.gfMap(lst)))'==2 ); %capillary nodes
lstNa = find( segVesType(im2.nodeSegN(im2.VesFlux.gfMap(lst)))'~=2 & segVesType(im2.nodeSegN(im2.VesFlux.gfMap(lst)))'~=3); %artery nodes
boo = ismember(im2.Mesh.boundary(:,1:3),lst(lstNav));
lstMav = find(sum(boo,2)==3); %arteries and veins mesh boundaries
boo = ismember(im2.Mesh.boundary(:,1:3),lst(lstNc));
lstMc = find(sum(boo,2)==3); %capillary mesh boundaries
boo = ismember(im2.Mesh.boundary(:,1:3),lst(lstNa));
lstMa = find(sum(boo,2)==3); %arteries mesh boundaries
lstAll=find(im2.Mesh.boundary(:,end)==0); %all nodes (a, c, v)


%%%% Display everything %%%%%%%%

if 0
    f10=figure(10);
    dim_fig=get(f10,'Position');
    dim_fig(3)=2.0*dim_fig(3);
    dim_fig(4)=1.0*dim_fig(4);
    set(f10,'Position',dim_fig);
    cm=colormap;
    subplot(1,2,1)
    for ii=4:size(segDilationDyn,2)
        plot(tFlow,segDilationDyn(:,ii),'color',cm(floor(64*ii/size(segDilationDyn,2)),:),'linewidth',2);hold on;
    end
    xlabel('time (s)','fontsize',18)
    title('Arterial dilation (%)','fontsize',18)
    axis square;
    legend('pial surface','diving 0-100um','diving 100-200um','diving 200-300um','diving 300-400um','diving 400-500um','diving 500-600um','diving >600um','1st br <150um','1st br >150um','2,3,4th br <150um','2,3,4th br >150um')
    subplot(1,2,2)
    trisurf( im2.Mesh.boundary(lstMa,1:3), im2.Mesh.node(:,1), im2.Mesh.node(:,2), -im2.Mesh.node(:,3), foo, 'linestyle','none','facecolor','flat');
    axis image
    view(22,34)
    colorbar
end

%%%% Run Flow-Volume dynamics %%%%%%%%

if 1
    
    %list of inflowing and outflowing segments
    endNodeIdx=find(im2.nB==1);    
    im2.nodeBCType(endNodeIdx)=1; %constant pressure on all end nodes
    im2.nodeBC=im2.nodePressure; %use pressure from baseline flow calculation results
    
    
    climvalueFlow=[-50 50];
    climvalueVol=[-10 10];
    climvaluePressure=[20 70];
    climvalueDilation=[-2 15];
    print_flag=8000; %number represent #of frame skipped between display
    [Ftime,Ptime,Vtime,Fedges,segNodeMap]=flowCircuitEq_dyn( im2, tFlow, seg_dvdt_t, segDilationDynType, lstAll, climvalueVol, climvalueFlow, climvalueDilation, climvaluePressure, lstSegActive, lstSegPassive, print_flag );
    save([savefolder,'dilate_vessel_NCES_Sigmoid',num2str(VesselDiPercent),'.mat'],'Ftime','Ptime','Vtime','Fedges','tFlow');
    
end






