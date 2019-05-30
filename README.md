# VAN_modeling

The vascular anatomical network (VAN) model computes BOLD signals of mouse microvascular stacks obtained from in-vivo two photon microscopy measurements. With this bottom-up model, you can explore different physiological effects on the oxygen distribution, the blood flow distribution and BOLD signals.

Here we provide two versions of instructions for users. The ‘Short Introduction to VAN’ outlined the steps to call the functions, which should be sufficient for people who are only interested in obtaining the results. The ‘Not So Short Introduction to VAN’ includes a review of VAN modeling and more details, which is for researchers who need more information or want to integrate our data with other computational methods. 

Any comments are welcome, and we would also love to discuss potential research topics. Emails with any questions can be sent to Xiaojun Cheng xcheng17@bu.edu, David Boas dboas@bu.edu. Most importantly, have fun with VAN!




 
Figure 1. Three-dimensional rendering of the six vascular stacks acquired with two-photon microscopy. 

## Citations
Cheng, X., Berman, A.J., Polimeni, J.R., Buxton, R.B., Gagnon, L., Devor, A., Sakadžić, S. and Boas, D.A., 2019. Dependence of the MR signal on the magnetic susceptibility of blood studied with models based on real microvascular networks. Magnetic resonance in medicine.

Gagnon, L., Sakadžić, S., Lesage, F., Pouliot, P., Dale, A.M., Devor, A., Buxton, R.B., Boas, D.A., 2016. Validation and optimization of hypercapnic-calibrated fMRI from oxygen-sensitive two-photon microscopy. Phil Trans R Soc B 371, 20150359.

Gagnon, L., Sakadžić, S., Lesage, F., Mandeville, E.T., Fang, Q., Yaseen, M.A., Boas, D.A., 2015a. Multimodal reconstruction of microvascular-flow distributions using combined two-photon microscopy and Doppler optical coherence tomography. Neurophotonics 2, 015008. 

Gagnon, L., Sakadžić, S., Lesage, F., Musacchia, J.J., Lefebvre, J., Fang, Q., Yücel, M.A., Evans, K.C., Mandeville, E.T., Cohen-Adad, J., Polimeni, J.R., Yaseen, M.A., Lo, E.H., Greve, D.N., Buxton, R.B., Dale, A.M., Devor, A., Boas, D.A., 2015. Quantifying the microvascular origin of BOLD-fMRI from first principles with two-photon microscopy and an oxygen-sensitive nanoprobe. J. Neurosci. Off. J. Soc. Neurosci. 35, 3663–3675.

Fang, Q., Sakadžić, S., Ruvinskaya, L., Devor, A., Dale, A.M., Boas, D.A., 2008. Oxygen advection and diffusion in a three dimensional vascular anatomical network. Opt. Express 16, 17530.






Short Introduction to VAN

	Download the zip file. Data of 6 VANs and scripts are included in the VANmodel_0.0.zip file
	Unzip all the files to a folder and copy the folder name. For example, in the below case, the folder name is ‘C:\Users\xcheng17\Documents\VAN\’
 
	Open the Matlab script VANmodelrun.m. Make sure you change ‘myfolder1’ to the folder you just saved your files. 
	If you are using Windows instead of Linux, make sure you change all the ‘/’ in all the scripts to ‘\’;
	You can specify the parameters you want, such as the magnetic field strength B0.
	Run the script. All the results, including intermediate results, are saved in the ‘data/mouse#/results/’ folder. Here the number ‘#’ is a number from 1 to 6 as you specified in the variable ‘mouseindex’.
	Do something else. The computation takes about 12-24 hours, depending on your computer speed.
	The gradient echo BOLD signal is the ‘BOLD_Gx’ variable in ‘MCBOLD_##.mat’.




Not So Short Introduction to VAN
Data
There are six folders in /data/ folder, each with vasculature graphs for a single mouse. 
Each mouse folder contains graph.seed, brOrder.mat, and Parameters.mat. These are the inputs for the model. The computed results from each module are saved in the /mouse#/results/ folder. 
The graph.seed file contains the graph of the vasculature. The information is contained in the im2 structure. Here’s a list of the most important fields of the im2 structure. To load the file into Matlab, go to the folder of the file and simply type ‘load('graph.seed', '-mat')’. If you are only interested in the vasculature, this ‘im2’ structure is what you are looking for.

im2.nodePos: the position of nodes. E.g. im2.nodePos(1,:) gives the x, y, z coordinates of the first node in μm.
im2.nodeEdges: the edges connecting the nodes. E.g. for the first mouse im2.nodeEdges(1,:)=15 41, which means that the 1st edge connects the 15th node and 41st node.
im2.nodeDiam: the diameter of the node in μm.
im2.nodeType: 1: arterioles; 2: capillaries; 3: venules.


Code
Mesh
The VAN_setup_mesh.m in the /code/ folder generates mesh, using iso2mesh toolbox. The input file is graph.seed and the output file is mesh.mat in /mouse#/results/ folder.
Vessel Dilation 
Input: mesh.mat, brOrder.mat.
Output: dilate_vessel_NCES-Sigmoid#.mat.
This module calculates flow dynamics with arterial dilation. The code used is VAN_calFlow_sigmoid.m.
Advection-diffusion solver
Input: mesh.mat, baseline CMRO2.
Output: advection_CMRO2_time_dilation ## ms.mat
The code to use is VAN_advection_setup.m. 
The advection-diffusion module is the most time-consuming part. Intermediate results are saved in /mouse#/results/. 
For example, if you have set the ‘duration_simulation’ time to be 18, but due to computer crash or some other reasons, you only see 16 files as below
 
You can use the function VAN_advection_setup_breaktoparts to complete the computation instead of rerun the simulation from the very beginning.
 


Monte-Carlo simulations
Input includes all the computed files from previous steps and the acquisition parameters of magnetic field strength and its orientation B0, gradient field Gx, echo time TE. 
The final result is BOLD_GE, which is the BOLD signal. 

Review of VAN modeling
Review of vascular anatomical network (VAN) modeling 
The VAN model is a bottom-up model that computes oxygen distribution, blood flow and BOLD signal within a micro-vascular network. 
Raw data obtained from experiments are converted to a graph which is a mathematical object of nodes and segments representing the vascular network. We have used a suite of custom-designed tools in MATLAB, with manual corrections applied until all segments in the field of view become interconnected. Vessel size was estimated at each graph node by thresholding the image at a low value of ~2%. Vessel types and branching orders are manually obtained by following them from the pial surface. Pial arterioles and venules are identified based on PO2 measurements and their morphology. 
Mesh generation is the first step of VAN modeling, providing meshes for the finite element method used to compute oxygen distribution. The toolbox we use is iso2mesh (Fang and Boas, 2009).
Oxygen distribution is obtained by solving the advection-diffusion equation 
(∂C_T)/∂t=v ⃗∙∇C_F-v ⃗∙∇C_B+∇∙(D_O2 ∇C_F )-OC,
where C_F and C_B are the free and bounded oxygen concentrations respectively; C_T=C_B+C_F denotes the total oxygen concentration, v ⃗ is the velocity, D_O2 is the oxygen diffusion coefficient and OC is the tissue oxygen consumption rate. Another equation that governs the equilibrium between C_F and C_B is
C_B=4C_Hb Hct〖SO〗_2 (C_F ),
where Hct is hematocrit; C_Hb is the hemoglobin concentration within a red blood cell, 〖SO〗_2 (C_F ) is the hemoglobin oxygen saturation which is in equilibrium with the concentration of free oxygen. This equation is solved by our hybrid model including a graph-based 1D oxygen advection within the vessel, a 1D oxygen flux conservation model across the vessel wall, and an finite-element-based 3D oxygen diffusion model within the tissue (Fang et al., 2008).
The VAN model computes blood flow within the vasculature from the flow-pressure relationship
P_1-P_2=F_12 R_12,
where P_1-P_2  is the pressure drop of the segment, F_12 is the blood flow and R_12 is the resistance which is inversely proportional to the square of the volume 1/〖V_12〗^2 (Boas et al., 2008). When the volume V_12 changes, pressure and blood flow distributions are recomputed to modulate neural activation with arterial dilation.
The Monte-Carlo method is used to compute the BOLD signal. The volume susceptibility shift is
∆χ=〖δχ〗_0 Hct(1-〖SO〗_2 ),
where  〖δχ〗_0 is the susceptibility difference between fully oxygenated and fully deoxygenated haemoglobin. The value of the haematocrit Hct is 0.3 in capillaries and 0.4 in arteries and veins. The magnetic field inhomogeneity ΔB(x ⃗) is computed by convolving ∆χ with the geometrical factor for the magnetic field induced by a unit cube
〖ΔB〗_cube=((2/π) a^3)/r^3 (3〖cos〗^2 θ-1)B_0.
Besides the magnetic field inhomogeneity, intrinsic T_2^* variations are also included in the modeling. It is very hard to model intravascular signals from microscopic details. We estimate T_2^* values inside vessels from experimental measurements (Uludağ et al., 2009)
T_(2,vessel)^* 〖=(A+〖C(1-〖SO〗_2 )〗^2)〗^(-1),
where A and C are constants that depend on the external magnetic field B_0 as in Table 1.
Table S1. Constants for T_2^* within vessels
Magnetic Field strength		A		C
B_0≤1.5T		6.5		25
〖1.5T<B〗_0≤3T		13.8		181
〖3T<B〗_0≤4T		30.4		262
〖4T<B〗_0≤4.7T		41		319
〖4.7T<B〗_0		100		500
In the tissue, T_2 and T_2^* are computed as
T_(2,tissue) 〖=(1.74∙B_0+7.77)〗^(-1),
T_(2,tissue)^* 〖=(3.74∙B_0+9.77)〗^(-1).
The diffusion of protons is modeled with a Monte-Carlo method. The initial positions of 〖10〗^7 protons are randomly distributed within the 3D volume. The diffusion coefficient of the protons is D=1×〖10〗^(-5)  〖cm〗^2 s^(-1), and the time step dt=0.2 ms. The position of each proton (x_1,x_2,x_3) after time dt is
x_i^'=x_i+N(0,2Ddt),i=1,2,3,
where N(0,2Ddt) is the normal distribution. Protons cannot move across the vessel membrane in our simulations. The phase increment for a proton at each dt is
〖δϕ〗_intra=-T_(2,vessel)^* dt/i ,
〖δϕ〗_intra=-T_(2,tissue)^* dt/i +γΔB(x ⃑)dt .
Here γ=2.675*〖10〗^5 rad/T/ms is the hydrogen proton precession frequency. The MR signal is computed as
S(t)=Re{1/N∑_1^N▒e^(〖iϕ〗_n (t)) }.
For simulation of the gradient echo MR signal, a gradient magnetic field ΔG_x x is turned on at TE/3 and flipped at 2TE/3. For the spin echo signal, the imaginary part of the phase is inverted at TE/2, ϕ_n (TE/2)=conj(ϕ_n (TE/2)). Finally, the BOLD signal is obtained as
BOLD=S(TE).
We can compute the BOLD signal for various oxygen consumption rates and blood flow distributions to model activation.

Notes
	VANmodelrun.m contains information to call different modules above. 
	Matlab versions beyond 2017a currently cannot calculate a large VAN. This is related to the algorithm of ‘\’ operator. We will work on this problem in the future.
	The real 〖CMRO〗_2 increase can be obtained from the increase of OCoutput variable in the advection_##_##_##ms.mat file from the output of the advection-diffusion solver.  This value can be different from the input 〖CMRO〗_2 change, which is currently a discrepency of our model. This difference is large for mouse2. This mouse was the first animal we dealth with and graphing may not be perfect. So if you are doing activation, you can avoid using this animal for now.


References
Boas, D.A., Jones, S.R., Devor, A., Huppert, T.J., Dale, A.M., 2008. A vascular anatomical network model of the spatio-temporal response to brain activation. Neuroimage 40, 1116–1129.
Fang, Q., Boas, D.A., 2009. Tetrahedral mesh generation from volumetric binary and grayscale images, in: 2009 IEEE International Symposium on Biomedical Imaging: From Nano to Macro. Presented at the 2009 IEEE International Symposium on Biomedical Imaging: From Nano to Macro, pp. 1142–1145. doi:10.1109/ISBI.2009.5193259
Fang, Q., Sakadžić, S., Ruvinskaya, L., Devor, A., Dale, A.M., Boas, D.A., 2008. Oxygen advection and diffusion in a three dimensional vascular anatomical network. Opt. Express 16, 17530.
Uludağ, K., Müller-Bierl, B., Uğurbil, K., 2009. An integrative model for neuronal activity-induced signal changes for gradient and spin echo functional imaging. NeuroImage 48, 150–165. doi:10.1016/j.neuroimage.2009.05.051


