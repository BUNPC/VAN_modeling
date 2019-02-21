
% This is a prototype for MCproton that will go from a pO2 graph
% distribution to a BOLD time course i.e. generation of the delta_B volume
% and MC simulation of protons.
%
% 11/21/2012 by L. Gagnon (adapted from F. Lesage's code)
% 
% 01/09/2013: - added volume changes when building SO2 and T2 masks L. Gagnon
%
% 02/12/2013: - added B0 dependence on T2_tissu from Uludag 2009 fit
%             - made vessel almost impermeable to proton
%             - separated IV from EV signal
%
% 03/14/2013: -added B0 orientation as an input parameter
%
% 04/30/2013: -assume no volume change for simplicity of IV and EV
%
% 11/01/2015: - got rid of compartment specific stuff to save memory
%             - vectorized the construction of the mask, got rid of loops (from P. Pouliot's code)
%             - vectorized construction of pert_B

function MCproton_Hybrid8(B0,phi_angle,omega_angle,TE,Gx,ROIx,ROIy,ROIz,file_name,iFrame,volume_and_mesh_file,vol_iFrame,ref_file,ref_iFrame,vol_ref_iFrame,VesselDiPercent,savefolder,cmro2)



%convert string input to numbers (when necessary)
B0=str2num(B0);
phi_angle=str2num(phi_angle);
omega_angle=str2num(omega_angle);
TE=str2num(TE);
Gx=str2num(Gx);
% ROIx1=str2num(ROIx1);
% ROIx2=str2num(ROIx2);
% ROIy1=str2num(ROIy1);
% ROIy2=str2num(ROIy2);
% ROIz1=str2num(ROIz1);
% ROIz2=str2num(ROIz2);

iFrame=str2num(iFrame);
vol_iFrame=str2num(vol_iFrame);
ref_iFrame=str2num(ref_iFrame);
vol_ref_iFrame=str2num(vol_ref_iFrame);

%%%%% MC parameters %%%%%%%%
nprotons=1e7;
dt=0.2;     %msec
ntime_step=ceil(TE/dt);
% ROIx=[ROIx1 ROIx2];    %um
% ROIy=[ROIy1 ROIy2];    %um
% ROIz=[ROIz1 ROIz2];    %um

%Volume size and voxel size (anisotropic)
Hvox=1; %in um (isotropic because the Fourier kernel assumes a sphere!)


% %%%%%% Some constants %%%%%%%%%%%%%%%%
% 
%B0=7; % main Field (Tesla) now as an input parameter
dChi0 = 4*pi*0.264e-6; % suceptibility of deoxygenated blood (from  Christen et al 2011)   
gamma=2.675e5; % rad/Tesla/msec
%Gx=20; % Gradient for readout in mT/m
Gxum=Gx*1e-6*1e-3; % Now in T/um
Dcoeff=1; %Proton (water) diffusion coefficient(um2/msec)
%TE=30; %Echo time (ms)
a=Hvox(1)/2; %radius of the sphere contained in a voxel (that's why we need isotropic voxels)
Dvox=Dcoeff/(Hvox(1).^2);

%Hematocrit
alpha = 1.27e-15;  % Bunsen solubility coefficient umol / um^3 / Torr
Hct_A = 0.44; %arteries (from Griffeth et al 2011)
Hct_C = 0.33; %capillaries
Hct_V = 0.44; %veins
Hct = [Hct_A Hct_C Hct_V];

%T2 constants
T2_tissu=1000.*(1.74*B0+7.77)^(-1); %in msec (from Uludag 2009)
T2star_tissu=1000.*(3.74*B0+9.77)^(-1); %in msec (from Uludag 2009)

%T2_blood = 1000.*(1./( (2.74.*B0-0.6)+(12.67.*B0^2.*(1-SatO2).^2) )); %in msec (from Uludag 2009)



%%%%%%% load everything %%%%%%%%%%%%%%%


%load the mesh and volume changes time course
load(volume_and_mesh_file);
load([savefolder,'dilate_vessel_NCES_Sigmoid',num2str(VesselDiPercent),'.mat']);
im2.Vtime=Vtime;

if ~isfield(im2,'species')
    species=1; %set to human pO2-sO2 curve otherwise specified
else
    species=im2.species;
end
if ~exist('Vtime')
    if isfield(im2,'Vtime')
        Vtime=im2.Vtime;
    else
        
        error('Don''t find Vtime variable')
    end
end
RelVol = Vtime(:,vol_iFrame) ./ Vtime(:,1);
RelVol_ref = Vtime(:,vol_ref_iFrame) ./ Vtime(:,1);


%load pO2 baseline file
load(ref_file);


cg1=cg(:,ref_iFrame);

clear c cg cbg

%load pO2 file
load(file_name);

cg2=cg(:,iFrame);

clear c cg cbg

%mask dimension
nx=ceil(im2.nX.*im2.Hvox(1));
ny=ceil(im2.nY.*im2.Hvox(2));
nz=ceil(im2.nZ.*im2.Hvox(3));

%generate names for files to be saved
[file_path,save_name_str,file_ext] = fileparts(file_name);
% save_name=['MCproton_Hybrid8_' save_name_str sprintf('_B0_%2.1fT_ori_[%d_%d]_TE%dms_Gx%d_ROI_[%d %d %d %d %d %d]_frame_%d_volFrame_%d',B0,phi_angle,omega_angle,TE,Gx,ROIx1,ROIx2,ROIy1,ROIy2,ROIz1,ROIz2,iFrame,vol_iFrame) '.mat'];
% save_name_mask=['mask' save_name];


%%%%%%% Generate SatO2 and T2 volumes %%%%%%%%%%%%%%%

% pO2 values are on graph nodes so I will generate a mask with a given voxel
% specifications.


%masking parameters (this allows a mask with no gap)
maskSphereSpacing = 0.05; %in um
maxR = 90; %in um

nEdges = size(im2.nodeEdges,1);
nNodes = size(im2.nodePos_um,1);
nodePos_um = im2.nodePos_um;

SatO2mask_ini = zeros(ny,nx,nz);
SatO2mask = zeros(ny,nx,nz);

Hctmask_ini = zeros(ny,nx,nz);
Hctmask = zeros(ny,nx,nz);
Vm_base = zeros(ny,nx,nz);
Vm = zeros(ny,nx,nz);

%Here we go
%hWait = waitbar( 0, sprintf('Masking vessels...'));

nodeMasked = zeros(nNodes,1);
nodeMasked_ref = zeros(nNodes,1);

lstEdges = 1:nEdges;
tic



for iiE=1:length(lstEdges)
    
    if(mod(iiE,floor(length(lstEdges)/10))==0)
        display(sprintf('Mask %d %%completed',round(100*iiE/length(lstEdges))))
    end
    iE = lstEdges(iiE);
    %waitbar(iiE/length(lstEdges),hWait);
    i1 = im2.nodeEdges(iE,1);
    i2 = im2.nodeEdges(iE,2);

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % actual volume size and O2 sat
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if (nodeMasked(i1)==0 || nodeMasked(i2)==0)
        if ~im2.maskUseSegDiam
            r1 = im2.nodeDiam(i1)/2.*RelVol(im2.nodeSegN(i1)).^(1/2);
            r1 = min(max(r1,1),maxR);
            r2 = im2.nodeDiam(i2)/2.*RelVol(im2.nodeSegN(i2)).^(1/2);
            r2 = min(max(r2,1),maxR);
            
        else
            foo = max( im2.segDiam(im2.nodeSegN(i1)).*RelVol(im2.nodeSegN(i1)).^(1/2) , im2.segDiam(im2.nodeSegN(i2)).*RelVol(im2.nodeSegN(i2)).^(1/2) )/2;
            foo = min(max(foo,1),maxR);
            r1 = foo;
            r2 = foo;
            
        end
        p1 = round(nodePos_um(i1,:)./Hvox); %in voxels (Hvox in the specified voxel size, not the native voxel size)
        p2 = round(nodePos_um(i2,:)./Hvox);
        d12 = norm(p2-p1);
        dxyz = (p2-p1);%/d12;
        rd = (r2-r1);%/d12;
        nSteps = 1;
        stepLen = max(r1*maskSphereSpacing,1);
        if stepLen<d12
            dxyz = (p2-p1)*stepLen/d12;
            rd = (r2-r1)*stepLen/d12;
            nSteps = floor(d12/stepLen)+1;
        end

        lst = find(sum((nodePos_um-ones(nNodes,1)*p1).^2,2).^0.5<(r1*maskSphereSpacing) );
        nodeMasked(lst) = 1;
        lst = find(sum((nodePos_um-ones(nNodes,1)*p2).^2,2).^0.5<(r2*maskSphereSpacing) );
        nodeMasked(lst) = 1;

        p = p1; %in voxels
        r = r1; %in um
        flag = 1;
        while flag
            if norm(round(p)-p2)==0
                flag = 0;
            end
            pr = round(p);
            rTmp = min(r,maxR);
            rx = ceil(rTmp/Hvox); %in voxels
            ry = ceil(rTmp/Hvox);
            rz = ceil(rTmp/Hvox);
            if im2.maskCircleFlag
                rz = 0;
            end
            
            if d12~=0
                
                %the satO2 value is constant for each voxel so compute it only
                %one time
                SatO2_value = norm(p-p1)/d12.*so2_func(cg2(i2)./alpha,species) + norm(p-p2)/d12.*so2_func(cg2(i1)./alpha,species);%weighted average

                %Philippe's code
                [Y1,Y2,Y3] = ndgrid(-rx:rx,-ry:ry,-rz:rz);
                Y1r = Y1(:);
                Y2r = Y2(:);
                Y3r = Y3(:);
                Y = [Y1r Y2r Y3r];
                Y_norm = sqrt(sum(Y.^2,2));
                Y_norm_small = find(Y_norm<=r);
                iix = min(max(pr(1)+Y1r(Y_norm_small),1),nx);
                iiy = min(max(pr(2)+Y2r(Y_norm_small),1),ny);
                iiz = min(max(pr(3)+Y3r(Y_norm_small),1),nz);
                SatO2mask(sub2ind([ny,nx,nz],iiy, iix, iiz)) = SatO2_value; % sub2ind, convert to a linear index
                %Hctmask(sub2ind([ny,nx,nz],iiy, iix, iiz)) = Hct(im2.segVesType(im2.edgeSegN(iiE)));
                Hctmask(sub2ind([ny,nx,nz],iiy, iix, iiz)) = Hct(ceil(im2.segVesType(im2.edgeSegN(iiE))));
                Vm(sub2ind([ny,nx,nz],iiy, iix, iiz)) = 1;
                
%                 for iX = -rx:+rx
%                     for iY = -ry:+ry
%                         for iZ = -rz:+rz
%                             if norm([iX*Hvox iY*Hvox iZ*Hvox])<=r %all in um
%                                 iix = min(max(pr(1)+iX,1),nx);
%                                 iiy = min(max(pr(2)+iY,1),ny);
%                                 iiz = min(max(pr(3)+iZ,1),nz);
%                                 Vm(iiy,iix,iiz) = 1;
%                                 SatO2mask(iiy,iix,iiz) = SatO2_value;
%                                 Hctmask(iiy,iix,iiz) = Hct(im2.segVesType(im2.edgeSegN(iiE)));
%                                 
%                             end
%                         end
%                     end
%                 end
            
            
            end
            
            if flag
                p = p + dxyz;
                r = r + rd;
            end
            nSteps = nSteps - 1;
            if nSteps==0
                flag = 0;
            end
            if norm(round(p)-p2)==0
                flag = 0;
            end

        end
    end % end of check if node type updated
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % reference volume size and O2 sat
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if (nodeMasked_ref(i1)==0 || nodeMasked_ref(i2)==0)
        if ~im2.maskUseSegDiam
            r1 = im2.nodeDiam(i1)/2.*RelVol_ref(im2.nodeSegN(i1)).^(1/2);
            r1 = min(max(r1,1),maxR);
            r2 = im2.nodeDiam(i2)/2.*RelVol_ref(im2.nodeSegN(i2)).^(1/2);
            r2 = min(max(r2,1),maxR);
            
        else
            foo = max( im2.segDiam(im2.nodeSegN(i1)).*RelVol_ref(im2.nodeSegN(i1)).^(1/2) , im2.segDiam(im2.nodeSegN(i2)).*RelVol_ref(im2.nodeSegN(i2)).^(1/2) )/2;
            foo = min(max(foo,1),maxR);
            r1 = foo;
            r2 = foo;
            
        end
        p1 = round(nodePos_um(i1,:)./Hvox); %in voxels (Hvox in the specified voxel size, not the native voxel size)
        p2 = round(nodePos_um(i2,:)./Hvox);
        d12 = norm(p2-p1);
        dxyz = (p2-p1);%/d12;
        rd = (r2-r1);%/d12;
        nSteps = 1;
        stepLen = max(r1*maskSphereSpacing,1);
        if stepLen<d12
            dxyz = (p2-p1)*stepLen/d12;
            rd = (r2-r1)*stepLen/d12;
            nSteps = floor(d12/stepLen)+1;
        end

        lst = find(sum((nodePos_um-ones(nNodes,1)*p1).^2,2).^0.5<(r1*maskSphereSpacing) );
        nodeMasked_ref(lst) = 1;
        lst = find(sum((nodePos_um-ones(nNodes,1)*p2).^2,2).^0.5<(r2*maskSphereSpacing) );
        nodeMasked_ref(lst) = 1;

        p = p1; %in voxels
        r = r1; %in um
        flag = 1;
        while flag
            if norm(round(p)-p2)==0
                flag = 0;
            end
            pr = round(p);
            rTmp = min(r,maxR);
            rx = ceil(rTmp/Hvox); %in voxels
            ry = ceil(rTmp/Hvox);
            rz = ceil(rTmp/Hvox);
            if im2.maskCircleFlag
                rz = 0;
            end
            
            if d12~=0
                
                %the satO2 value is constant for each voxel so compute it only
                %one time
                SatO2_value_ref = norm(p-p1)/d12.*so2_func(cg1(i2)./alpha,species) + norm(p-p2)/d12.*so2_func(cg1(i1)./alpha,species);%weighted average
                
                %Philippe's code
                [Y1,Y2,Y3] = ndgrid(-rx:rx,-ry:ry,-rz:rz);
                Y1r = Y1(:);
                Y2r = Y2(:);
                Y3r = Y3(:);
                Y = [Y1r Y2r Y3r];
                Y_norm = sqrt(sum(Y.^2,2));
                Y_norm_small = find(Y_norm<=r);
                iix = min(max(pr(1)+Y1r(Y_norm_small),1),nx);
                iiy = min(max(pr(2)+Y2r(Y_norm_small),1),ny);
                iiz = min(max(pr(3)+Y3r(Y_norm_small),1),nz);
                SatO2mask_ini(sub2ind([ny,nx,nz],iiy, iix, iiz)) = SatO2_value_ref;
                %Hctmask_ini(sub2ind([ny,nx,nz],iiy, iix, iiz)) = Hct(im2.segVesType(im2.edgeSegN(iiE)));
                
                Hctmask_ini(sub2ind([ny,nx,nz],iiy, iix, iiz)) = Hct(ceil(im2.segVesType(im2.edgeSegN(iiE))));
                
                Vm_base(sub2ind([ny,nx,nz],iiy, iix, iiz)) = 1;
                 
                
%                 for iX = -rx:+rx
%                     for iY = -ry:+ry
%                         for iZ = -rz:+rz
%                             if norm([iX*Hvox iY*Hvox iZ*Hvox])<=r %all in um
%                                 iix = min(max(pr(1)+iX,1),nx);
%                                 iiy = min(max(pr(2)+iY,1),ny);
%                                 iiz = min(max(pr(3)+iZ,1),nz);
%                                 SatO2mask_ini(iiy,iix,iiz) = SatO2_value_ref;     
%                                 Hctmask_ini(iiy,iix,iiz) = Hct(im2.segVesType(im2.edgeSegN(iiE)));
%                                 Vm_base(iiy,iix,iiz) = 1;
% 
%                             end
%                         end
%                     end
%                 end
            
            
            end
            
            if flag
                p = p + dxyz;
                r = r + rd;
            end
            nSteps = nSteps - 1;
            if nSteps==0
                flag = 0;
            end
            if norm(round(p)-p2)==0
                flag = 0;
            end

        end
    end % end of check if node type updated
    
    
end
%close(hWait)
toc



%T2 volume (in msec)
T2mask_ini = T2_tissu.*ones(ny,nx,nz);
T2mask_ini(Hctmask_ini~=0) = 1000.*(1./( (2.74.*B0-0.6)+(12.67.*B0^2.*(1-SatO2mask_ini(Hctmask_ini~=0)).^2) ));
T2mask = T2_tissu.*ones(ny,nx,nz);
T2mask(Hctmask~=0) = 1000.*(1./( (2.74.*B0-0.6)+(12.67.*B0^2.*(1-SatO2mask(Hctmask~=0)).^2) ));

%T2* volume (in msec)
%intravascular relaxation rate
%from Zhao 2007 and Silverstein 2003
if B0<=1.5
    As=6.5;
    Cs=25;
elseif (B0>1.5 && B0<=3)
    As=13.8;
    Cs=181;
elseif (B0>3 && B0<=4)
    As=30.4;
    Cs=262;
elseif (B0>4 && B0<=4.7)
    As=41;
    Cs=319;
elseif B0>4.7
    As=100;
    Cs=500; %very large number so that decays faster than TE
end

T2starmask_ini = T2star_tissu.*ones(ny,nx,nz);
T2starmask_ini(Hctmask_ini~=0) = 1000.*(1./(As+Cs.*(1-SatO2mask_ini(Hctmask_ini~=0)).^2));
T2starmask = T2star_tissu.*ones(ny,nx,nz);
T2starmask(Hctmask~=0) = 1000.*(1./(As+Cs.*(1-SatO2mask(Hctmask~=0)).^2));


%dChi volume
dChimask_ini = dChi0.*Hctmask_ini.*(1-SatO2mask_ini);
dChimask = dChi0.*Hctmask.*(1-SatO2mask);
clear Hctmask Hctmask_ini

%crop volumes to keep only ROI (x and y are inverted)
xidx1=max(round(ROIy(1)),1);
xidx2=min(round(ROIy(2)),nx);
yidx1=max(round(ROIx(1)),1);
yidx2=min(round(ROIx(2)),ny);
zidx1=max(round(ROIz(1)),1);
zidx2=min(round(ROIz(2)),nz);
xidx2=floor((xidx2-xidx1+1)/2)*2+(xidx1-1);
yidx2=floor((yidx2-yidx1+1)/2)*2+(yidx1-1);
zidx2=floor((zidx2-zidx1+1)/2)*2+(zidx1-1);

dChimask_ini = dChimask_ini(xidx1:xidx2,yidx1:yidx2,zidx1:zidx2);
dChimask = dChimask(xidx1:xidx2,yidx1:yidx2,zidx1:zidx2);

T2starmask_ini = T2starmask_ini(xidx1:xidx2,yidx1:yidx2,zidx1:zidx2);
T2starmask = T2starmask(xidx1:xidx2,yidx1:yidx2,zidx1:zidx2);

T2mask_ini = T2mask_ini(xidx1:xidx2,yidx1:yidx2,zidx1:zidx2);
T2mask = T2mask(xidx1:xidx2,yidx1:yidx2,zidx1:zidx2);

Vm = Vm(xidx1:xidx2,yidx1:yidx2,zidx1:zidx2);
Vm_base = Vm_base(xidx1:xidx2,yidx1:yidx2,zidx1:zidx2);

[nx ny nz]=size(dChimask);



%%%%%%%%% Construct delta_B map %%%%%%%%%%%%%%
%(we assume cube is repeated periodically in space, no padding in FFT)

fft_succep_volume=fftshift(fftn(fftshift(dChimask)));
fft_succep_volume_ini=fftshift(fftn(fftshift(dChimask_ini)));
clear dChimasl dChimask_ini

origin=[nx/2-0.5,ny/2-0.5,nz/2-0.5].*Hvox(1); %in um
pert_B=zeros(nx,ny,nz);
r0 = [ cos(phi_angle*pi/180)*sin(omega_angle*pi/180) sin(phi_angle*pi/180)*sin(omega_angle*pi/180) cos(omega_angle*pi/180) ]; %unit vector oriented along the B-field
r0_norm=norm(r0); %should be one

%Philippe's code
%List of all vectors through the origin
[X1,X2,X3] = ndgrid(1:nx,1:ny,1:nz);
X1r = X1(:);
X2r = X2(:);
X3r = X3(:);
X = [X1r X2r X3r]-ones(size(X1r,1),1)*origin;
r_norm = sqrt(sum(X.^2,2));
costheta=dot(X,ones(size(X1r,1),1)*r0,2)./(r_norm.*r0_norm);
pert_B=B0*2/pi*a^3.*(r_norm).^(-3).*(3*costheta.^2-1);
pert_B(isinf(pert_B)) = 0; %Does this ever happen?
pert_B(isnan(pert_B)) = 0;
pert_B = reshape(pert_B,[nx,ny,nz]);

% for ii=1:size(pert_B,1)
%     for jj=1:size(pert_B,2)
%         for kk=1:size(pert_B,3)
%             r=[ii,jj,kk].*Hvox(1)-origin; %in um
%             r_norm=norm(r); %in um
%             if r_norm==0
%                 pert_B(ii,jj,kk)=0;
%             else
%                 costheta=dot(r,r0)./(r_norm.*r0_norm);
%                 pert_B(ii,jj,kk)=B0*2/pi*a^3/(r_norm)^3*(3*costheta^2-1);
%             end
%         end
%     end
% end


fft_pert_B=fftshift(fftn(fftshift(pert_B)));
clear pert_B
delta_B_ini=real(ifftshift(ifftn(ifftshift(fft_pert_B.*fft_succep_volume_ini))));
delta_B=real(ifftshift(ifftn(ifftshift(fft_pert_B.*fft_succep_volume))));
clear fft_succep_volume fft_succep_volume_ini

%%%%%%%%% Run proton Monte Carlo %%%%%%%%%%%%

phase_GE_base=zeros(nprotons,1);
phase_SE_base=zeros(nprotons,1);
signal_GE_base=zeros(1,ntime_step);
signal_SE_base=zeros(1,ntime_step);
signal_GE_IV_base=zeros(1,ntime_step);
signal_SE_IV_base=zeros(1,ntime_step);
signal_GE_EV_base=zeros(1,ntime_step);
signal_SE_EV_base=zeros(1,ntime_step);


phase_GE=zeros(nprotons,1);
phase_SE=zeros(nprotons,1);

signal_GE=zeros(1,ntime_step);
signal_SE=zeros(1,ntime_step);
signal_GE_IV=zeros(1,ntime_step);
signal_SE_IV=zeros(1,ntime_step);
signal_GE_EV=zeros(1,ntime_step);
signal_SE_EV=zeros(1,ntime_step);


% Generate initial position for protons inside the ROI
%protons_pos=floor(nx*rand(nprotons,3))+0.5;
protons_pos_x=floor(nx*rand(nprotons,1))+0.5;  %this is in microns
protons_pos_y=floor(ny*rand(nprotons,1))+0.5;  %this is in microns
protons_pos_z=floor(nz*rand(nprotons,1))+0.5;  %this is in microns
protons_pos=[protons_pos_x protons_pos_y protons_pos_z];
protons_pos_base=[protons_pos_x protons_pos_y protons_pos_z];    %need this because proton must stay inside vessel at baseline too


%Spin Echo params
half_echo_index=round(TE/(2*dt));
echo_index=min(round(TE/dt),ntime_step);
start1_se_index=round(TE/(4*dt)-TE/(8*dt));
stop1_se_index=round(TE/(4*dt)+TE/(8*dt));
start2_se_index=round(TE/(dt)-TE/(4*dt));
stop2_se_index=round(TE/dt+TE/(4*dt));

%Gradient Echo params
start_ge_index=round(TE/(3*dt));
flip_ge_index=round(2*TE/(3*dt));
stop_ge_index=round(4*TE/(3*dt));


ge_gradient_weight=0;
se_gradient_weight=0;
tic
for i=1:ntime_step
    if(mod(i,(ntime_step/10))==0)
        display(sprintf('MCproton %d %%completed',round(100*i/ntime_step)))
    end
    
    %Gradient in readout for GE
    if(start_ge_index==i)
        ge_gradient_weight=1;
    elseif (flip_ge_index==i)
        ge_gradient_weight=-1;
    elseif (stop_ge_index==i)
        ge_gradient_weight=0;
    end
    
    %Gradient in readout for SE
    if(start1_se_index==i)
        se_gradient_weight=1;
    elseif (stop1_se_index==i)
        se_gradient_weight=0;
    elseif (start2_se_index==i)
        se_gradient_weight=1;
    elseif (stop2_se_index==i)
        se_gradient_weight=0;
    end
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%
    % activated state
   
    %identify extravascular proton
    indices=ceil(protons_pos/Hvox(1));
    indices(:,1)=max(min(indices(:,1),nx),1);
    indices(:,2)=max(min(indices(:,2),ny),1);
    indices(:,3)=max(min(indices(:,3),nz),1);
    proton_Vm=Vm(sub2ind(size(Vm),indices(:,1),indices(:,2),indices(:,3))); 
    
    %move proton temporarly
    sigma=sqrt(2*Dvox*dt);
    rnum=sigma*randn(nprotons,3);
    protons_pos_new_temp=protons_pos+rnum;
    
    %(make sure IV stays IV and EV stays EV)
    indices=ceil(protons_pos_new_temp/Hvox(1));
    indices(:,1)=max(min(indices(:,1),nx),1);
    indices(:,2)=max(min(indices(:,2),ny),1);
    indices(:,3)=max(min(indices(:,3),nz),1);
    proton_Vm_new=Vm(sub2ind(size(Vm),indices(:,1),indices(:,2),indices(:,3))); 
    
    %find the ones who crossed vessel boundary (try 100 times otherwise let it go cross the boundary)
    crossing_list=find(proton_Vm~=proton_Vm_new);
    flag_cpt=0;
    while ( ~isempty(crossing_list) && flag_cpt<=100 )
        
        rnum(crossing_list,:)=sigma*randn(length(crossing_list),3);
        protons_pos_new_temp=protons_pos+rnum;
        
        indices=ceil(protons_pos_new_temp/Hvox(1));
        indices(:,1)=max(min(indices(:,1),nx),1);
        indices(:,2)=max(min(indices(:,2),ny),1);
        indices(:,3)=max(min(indices(:,3),nz),1);
        proton_Vm_new=Vm(sub2ind(size(Vm),indices(:,1),indices(:,2),indices(:,3))); 
        crossing_list=find(proton_Vm~=proton_Vm_new);
    
        flag_cpt=flag_cpt+1;
        
    end
    
    %move proton accordingly
    protons_pos=protons_pos+rnum;
    indices=ceil(protons_pos/Hvox(1));
    indices(:,1)=max(min(indices(:,1),nx),1);
    indices(:,2)=max(min(indices(:,2),ny),1);
    indices(:,3)=max(min(indices(:,3),nz),1);
    
    %update IV and EV lists
    proton_Vm=Vm(sub2ind(size(Vm),indices(:,1),indices(:,2),indices(:,3)));
    IV_list=find(proton_Vm==1);  % here intravascular and extravascular signal are seperated
    EV_list=find(proton_Vm~=1);
    %display(sprintf('%1.4f percent of protons were intravascular',100.*length(IV_list)./nprotons))
    
    
    %%%%%%%%%%%%%%%%%%%%
    % baseline state
   
    %identify extravascular proton
    indices_base=ceil(protons_pos_base/Hvox(1));
    indices_base(:,1)=max(min(indices_base(:,1),nx),1);
    indices_base(:,2)=max(min(indices_base(:,2),ny),1);
    indices_base(:,3)=max(min(indices_base(:,3),nz),1);
    proton_Vm_base=Vm_base(sub2ind(size(Vm_base),indices_base(:,1),indices_base(:,2),indices_base(:,3))); 
    
    %move proton temporarly
    sigma=sqrt(2*Dvox*dt);
    rnum_base=sigma*randn(nprotons,3);
    protons_pos_base_new_temp=protons_pos_base+rnum_base;
    
    %(make sure IV stays IV and EV stays EV)
    indices_base=ceil(protons_pos_base_new_temp/Hvox(1));
    indices_base(:,1)=max(min(indices_base(:,1),nx),1);
    indices_base(:,2)=max(min(indices_base(:,2),ny),1);
    indices_base(:,3)=max(min(indices_base(:,3),nz),1);
    proton_Vm_base_new=Vm_base(sub2ind(size(Vm_base),indices_base(:,1),indices_base(:,2),indices_base(:,3))); 
    
    %find the ones who crossed vessel boundary (try 100 times otherwise let it go cross the boundary)
    crossing_list=find(proton_Vm_base~=proton_Vm_base_new);
    flag_cpt=0;
    while ( ~isempty(crossing_list) && flag_cpt<=100 )
        
        rnum_base(crossing_list,:)=sigma*randn(length(crossing_list),3);
        protons_pos_base_new_temp=protons_pos_base+rnum_base;
        
        indices_base=ceil(protons_pos_base_new_temp/Hvox(1));
        indices_base(:,1)=max(min(indices_base(:,1),nx),1);
        indices_base(:,2)=max(min(indices_base(:,2),ny),1);
        indices_base(:,3)=max(min(indices_base(:,3),nz),1);
        proton_Vm_base_new=Vm_base(sub2ind(size(Vm_base),indices_base(:,1),indices_base(:,2),indices_base(:,3))); 
        crossing_list=find(proton_Vm_base~=proton_Vm_base_new);
    
        flag_cpt=flag_cpt+1;
        
    end
    
    %move proton accordingly
    protons_pos_base=protons_pos_base+rnum_base;
    indices_base=ceil(protons_pos_base/Hvox(1));
    indices_base(:,1)=max(min(indices_base(:,1),nx),1);
    indices_base(:,2)=max(min(indices_base(:,2),ny),1);
    indices_base(:,3)=max(min(indices_base(:,3),nz),1);
    
    %update IV and EV lists
    proton_Vm_base=Vm_base(sub2ind(size(Vm_base),indices_base(:,1),indices_base(:,2),indices_base(:,3)));
    IV_list_base=find(proton_Vm_base==1);
    EV_list_base=find(proton_Vm_base~=1);
    %display(sprintf('%1.4f percent of protons were intravascular at baseline',100.*length(IV_list_base)./nprotons))
    
    
    protons_B_base=delta_B_ini(sub2ind(size(delta_B_ini),indices_base(:,1),indices_base(:,2),indices_base(:,3)));
    protons_B=delta_B(sub2ind(size(delta_B),indices(:,1),indices(:,2),indices(:,3)));
    
    protons_T2_base=T2mask_ini(sub2ind(size(delta_B_ini),indices_base(:,1),indices_base(:,2),indices_base(:,3)));
    protons_T2=T2mask(sub2ind(size(delta_B),indices(:,1),indices(:,2),indices(:,3)));
    
    protons_T2star_base=T2starmask_ini(sub2ind(size(delta_B_ini),indices_base(:,1),indices_base(:,2),indices_base(:,3)));
    protons_T2star=T2starmask(sub2ind(size(delta_B),indices(:,1),indices(:,2),indices(:,3)));
    
    protons_gradient_GE=ge_gradient_weight*(protons_pos(:,1)-(nx/2)*Hvox(1))*Gxum;
    protons_gradient_SE=se_gradient_weight*(protons_pos(:,1)-(nx/2)*Hvox(1))*Gxum;
    protons_gradient_GE_base=ge_gradient_weight*(protons_pos_base(:,1)-(nx/2)*Hvox(1))*Gxum;
    protons_gradient_SE_base=se_gradient_weight*(protons_pos_base(:,1)-(nx/2)*Hvox(1))*Gxum;
    
    
    % Spin echo, 180 degrees phase shift
    if i==half_echo_index
      
        phase_SE_base=conj(phase_SE_base); 
        phase_SE=conj(phase_SE);
        
    end
    
    
    %%%%%%%%%%%%%%%
    % baseline signal (all vessels)
    
    %signal GE
    phase_GE_base(EV_list_base)=phase_GE_base(EV_list_base)+j*gamma*(protons_B_base(EV_list_base)+protons_gradient_GE_base(EV_list_base))*dt-(1./protons_T2_base(EV_list_base))*dt;
    phase_GE_base(IV_list_base)=phase_GE_base(IV_list_base)-(1./protons_T2star_base(IV_list_base))*dt;
    signal_GE_base(i)=abs(sum(exp(phase_GE_base))/nprotons);
    signal_GE_IV_base(i)=abs(sum(exp(phase_GE_base(IV_list_base)))/length(IV_list_base));
    signal_GE_EV_base(i)=abs(sum(exp(phase_GE_base(EV_list_base)))/length(EV_list_base));
   
    %signal SE
    phase_SE_base(EV_list_base)=phase_SE_base(EV_list_base)+j*gamma*(protons_B_base(EV_list_base)+protons_gradient_SE_base(EV_list_base))*dt-(1./protons_T2_base(EV_list_base))*dt;
    phase_SE_base(IV_list_base)=phase_SE_base(IV_list_base)-(1./protons_T2_base(IV_list_base))*dt;
    signal_SE_base(i)=abs(sum(exp(phase_SE_base))/nprotons);
    signal_SE_IV_base(i)=abs(sum(exp(phase_SE_base(IV_list_base)))/length(IV_list_base));
    signal_SE_EV_base(i)=abs(sum(exp(phase_SE_base(EV_list_base)))/length(EV_list_base));
    
    
    %%%%%%%%%%%%%%%
    % Activation signal (all vessels)
    
    %signal GE
    phase_GE(EV_list)=phase_GE(EV_list)+j*gamma*(protons_B(EV_list)+protons_gradient_GE(EV_list))*dt-(1./protons_T2(EV_list))*dt;
    phase_GE(IV_list)=phase_GE(IV_list)-(1./protons_T2star(IV_list))*dt;
    signal_GE(i)=abs(sum(exp(phase_GE))/nprotons);
    signal_GE_IV(i)=abs(sum(exp(phase_GE(IV_list)))/length(IV_list));
    signal_GE_EV(i)=abs(sum(exp(phase_GE(EV_list)))/length(EV_list));
      
    %signal SE
    phase_SE(EV_list)=phase_SE(EV_list)+j*gamma*(protons_B(EV_list)+protons_gradient_SE(EV_list))*dt-(1./protons_T2(EV_list))*dt;
    phase_SE(IV_list)=phase_SE(IV_list)-(1./protons_T2(IV_list))*dt;
    signal_SE(i)=abs(sum(exp(phase_SE))/nprotons);
    signal_SE_IV(i)=abs(sum(exp(phase_SE(IV_list)))/length(IV_list));
    signal_SE_EV(i)=abs(sum(exp(phase_SE(EV_list)))/length(EV_list));
    
    
    
    
end
toc

tMCproton=(1:ntime_step)*dt;
BOLD_GE_base=signal_GE_base(echo_index);
BOLD_SE_base=signal_SE_base(echo_index);
BOLD_GE_IV_base=signal_GE_IV_base(echo_index);
BOLD_SE_IV_base=signal_SE_IV_base(echo_index);
BOLD_GE_EV_base=signal_GE_EV_base(echo_index);
BOLD_SE_EV_base=signal_SE_EV_base(echo_index);

BOLD_GE=signal_GE(echo_index)
BOLD_SE=signal_SE(echo_index);
BOLD_GE_IV=signal_GE_IV(echo_index);
BOLD_SE_IV=signal_SE_IV(echo_index);
BOLD_GE_EV=signal_GE_EV(echo_index);
BOLD_SE_EV=signal_SE_EV(echo_index);

save([savefolder,'MCBOLD',num2str(VesselDiPercent),'_',num2str(cmro2),'.mat'],'tMCproton','TE','signal_GE','signal_SE','signal_GE_IV','signal_SE_IV','signal_GE_EV','signal_SE_EV',...
    'signal_GE_base','signal_SE_base','signal_GE_IV_base','signal_SE_IV_base','signal_GE_EV_base','signal_SE_EV_base',...
    'BOLD_GE','BOLD_SE','BOLD_GE_IV','BOLD_SE_IV','BOLD_GE_EV','BOLD_SE_EV','BOLD_GE_base','BOLD_SE_base','BOLD_GE_IV_base',...
    'BOLD_SE_IV_base','BOLD_GE_EV_base','BOLD_SE_EV_base') %this is a small file

%eval(sprintf('save %s tMCproton TE signal_GE signal_SE signal_GE_IV signal_SE_IV signal_GE_EV signal_SE_EV      signal_GE_base signal_SE_base signal_GE_IV_base signal_SE_IV_base signal_GE_EV_base signal_SE_EV_base   BOLD_GE BOLD_SE BOLD_GE_IV BOLD_SE_IV BOLD_GE_EV BOLD_SE_EV      BOLD_GE_base BOLD_SE_base BOLD_GE_IV_base BOLD_SE_IV_base BOLD_GE_EV_base BOLD_SE_EV_base',save_name)) %this is a small file
%eval(sprintf('save %s delta_B SatO2mask T2mask Hctmask',save_name_mask)) %this gives huge files! (~1GB) 

%save([savefolder,'MCBOLD',num2str(VesselDiPercent),'_',num2str(cmro2),'.mat'],'BOLD_GE_base','BOLD_GE') %this is a small file
end


