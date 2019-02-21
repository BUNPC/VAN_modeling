%%some advection results

function VAN_advection_setup_breaktoparts(CMRO2_baseline,duration_simulation,savefolder,VesselDiPercent,lastfileN)


% Load mesh
load([savefolder,'mesh.mat'])  ; 







im2.species=2;

im2.nIter = 200;
im2.nIterPerTrecord = 1;
im2.nAdvIterPerIter = 40;

im2.dt  = 5e-3;  % units of seconds
im2.Hct = 0.3;

im2.Do2 = 2.4e3;    % 2.4e-3 mm^2/s = 2.4e3 um^2/s

im2.hwall = 1;    % vessel wall thickness um
im2.Kves  = 2.23e-12; % Kves = 5e-8 uL mm^-1 s^-1 mmHg^-1 Popel1989
    % 1 L = 1000 cm^3 = 1e6 mm^3 = 1e15 um^3 = 1e6 uL
    % Kves = 5e-2 um^3 um^-1 s^-1 mmHg^-1
    % need to convert from vol of O2 to moles of O2
    % multiply by Pstandard / (R T)
    %     = 1 / (8.21e13 * 273)
    %     = 4.46e-17 mol / um^3
    % Pstandard = 1 atm
    % R = .0821 L atm K^-1 mol^-1
    %   = 8.21e13 um^3 atm K^-1 mol^-1
    % T = 273 K
    % Kves = 2.23e-18 mol um^-1 s^-1 mmHg^-1
    % Kves = 2.23e-12 umol um^-1 s^-1 mmHg^-1
im2.KvesArtFact=1;

% im2.po2a = 90;      % Torr at the arteriole side Input
im2.po2t = 10;      % initial conditions for po2 everywhere

%CMRO2 (trace shape is average of active dilation)
im2.OCo_umol_per_mL_per_min = CMRO2_baseline; % umol / cm^3 / min
im2.OCo = im2.OCo_umol_per_mL_per_min/1e12/60; % umol / um^3 / sec
tCMRO2 = 0;
im2.tOC = tCMRO2';

%This is for tissue or blood???, Michele found that it wasn't the same
im2.alpha = 1.27e-15;  % Bunsen solubility coefficient umol / um^3 / Torr

im2.Chb = 5.3e-12;            % 5.3e3 uM  red cell hemoglobin content (Beard2001)
% 5.3e3 umol/L = 5.3e-12 umol/um^3


%FLOW
tFlow = 0:0.025:18.1; % This 18.1 should not be hard coded. It needs to be longer than duration_simulation. Also, the 0.025 can be different than dt. It is purposefully more coarse.
im2.tFlow = tFlow'; %must be nT-by-1 at this point
im2.Fedges = im2.edgeFlow*ones(1,length(tFlow)); %must be nEdge-by-nT
im2.Ftime = [im2.segFlow; im2.segFlow]*ones(1,length(tFlow)); %must be 2*nSeg-by-nT
im2.segRelVol = ones(length(im2.segFlow),length(im2.tFlow)); %must be nSeg-by-nT

%display
nIterPerPlot = im2.nIter+1; %never display

if VesselDiPercent ~=0
% Load mesh
load([savefolder,'dilate_vessel_NCES_Sigmoid',num2str(VesselDiPercent),'.mat']);
tFlow=tFlow-tFlow(1);
im2.tFlow = tFlow';
im2.Fedges = Fedges;
im2.segRelVol = Vtime ./ (Vtime(:,1)*ones(1,length(tFlow)));
im2.Ftime = Ftime;

%add one point for the last itteration
dtFlow=mean(diff(tFlow));
im2.tFlow(end+1)=im2.tFlow(end)+dtFlow;
im2.Fedges(:,end+1)=im2.Fedges(:,end);
im2.segRelVol(:,end+1)=im2.segRelVol(:,end);
im2.Ftime(:,end+1)=im2.Ftime(:,end);
end

% % % BCs on inflowing nodes
% % Use branching order now. For artery use artery and for vein use vein. For
% % capillary, use vessel type corresponding to lowest branching order.
% % 5/11/2012
% 
% %                            Arteries          Veins
% %              br. order      pO2               pO2
% %
% TableFromSava = [ 0           102               45; %pial
%                   1            89               40;                                      
%                   2            78               36;
%                   3            72               33;
%                   4            70               32;
%                   5            58               30;
%                   6            47               30;
%                   7            38               30;
%                   8            35               30;
%                   9            35               30;
%                   10           35               30;
%                   ];
 
lstBC = find(im2.Adv.nodeIn==1);
nBC = length(lstBC);
BCvec=zeros(nBC,2);
for iBC=1:nBC
    BCvec(iBC,1)=lstBC(iBC);
    
    % set inflowing artery only and put the rest to 0
    if im2.nodeType(lstBC(iBC))==1 %artery
        BCvec(iBC,2)=102;
    else
        BCvec(iBC,2)=0; %otherwise zeros because flow is zero anyway!
    end
    
%     edgeIdx=find(im2.nodeEdges(:,1)==lstBC(iBC) | im2.nodeEdges(:,2)==lstBC(iBC));
%     brOrderArt=im2.edgeBRorderArt(edgeIdx);
%     brOrderVeins=im2.edgeBRorderVeins(edgeIdx);
%     
%     if im2.nodeType(lstBC(iBC))==1 %artery
%         BCvec(iBC,2)=TableFromSava(min(brOrderArt+1,11),2);
%     elseif im2.nodeType(lstBC(iBC))==3 %veins
%         BCvec(iBC,2)=TableFromSava(min(brOrderVeins+1,11),3);
%     elseif im2.nodeType(lstBC(iBC))==2 %capillary, take minimum br. Order
%         BCvec(iBC,2)=TableFromSava(min(min(brOrderArt,brOrderVeins)+1,11),3);
%     end
       
end

%######################################
% Run the simulation
%######################################
% if lastfileN~=0
%load([savefolder,'FluxMatrices.mat']);
   load([savefolder,'advection_',num2str(CMRO2_baseline),'_',num2str(lastfileN*1000),'_',num2str(VesselDiPercent),'ms.mat'])
     

for iSim=(lastfileN+1):duration_simulation
    
    
    tStart=iSim-1;
    if iSim==1
        [c,cg,cbg,t,AmatFlux,BmatFlux,OCoutput] = VAN_advection_run( im2, im2.nIter, nIterPerPlot, [], [], [], [], [], tStart, BCvec );
        save([savefolder,'FluxMatrices.mat'],'AmatFlux','BmatFlux');
    else
        c1=c(:,end);
        cg1=cg(:,end);
        cbg1=cbg(:,end);
        [c,cg,cbg,t,AmatFlux,BmatFlux,OCoutput] = VAN_advection_run( im2, im2.nIter, nIterPerPlot, AmatFlux, BmatFlux, c1, cg1, cbg1, tStart, BCvec );
    end
    
    %eval(sprintf('save advection_%1.1f_%dms.mat c cg cbg OCoutput',CMRO2_baseline,iSim*1000))
     save([savefolder,'advection_',num2str(CMRO2_baseline),'_',num2str(iSim*1000),'_',num2str(VesselDiPercent),'ms.mat'],'c','cg','cbg','OCoutput','t','AmatFlux','BmatFlux');
end

return