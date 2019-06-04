# VAN_modeling

The vascular anatomical network (VAN) model computes *oxygen distribution*, *blood flow distribution* and *BOLD fMRI signals* within mouse  microvascular stacks obtained in-vivo from two-photon microscopy.

Here we provide instructions to run the code and [descriptions of the model](https://github.com/BUNPC/VAN_modeling/tree/master/Model_Description).

Any comments are welcome, and we would also love to discuss potential research topics. Emails with any questions can be sent to Xiaojun Cheng xcheng17@bu.edu, David Boas dboas@bu.edu. Most importantly, have fun with VAN!

![Figure](VANs.png)


 
Figure 1. Three-dimensional rendering of the six vascular stacks acquired with two-photon microscopy. 

## Citations
Cheng, X., Berman, A.J.L.J., Polimeni, J.R., Buxton, R.B., Gagnon, L., Devor, A., Sakadžić, S., and Boas, D.A., “Dependence of the MR signal on the magnetic susceptibility of blood studied with models based on real microvascular networks.,” Magnetic resonance in medicine (2019).

Gagnon, L., Sakadžić, S., Lesage, F., Pouliot, P., Dale, A.M., Devor, A., Buxton, R.B., and Boas, D.A., “Validation and optimization of hypercapnic-calibrated fMRI from oxygen-sensitive two-photon microscopy.,” Philosophical transactions of the Royal Society of London. Series B, Biological sciences 371(1705), 20150359 (2016).

Gagnon, L., Sakadžić, S., Lesage, F., Musacchia, J.J., Lefebvre, J., Fang, Q., Yücel, M.A., Evans, K.C., Mandeville, E.T., et al., “Quantifying the microvascular origin of BOLD-fMRI from first principles with two-photon microscopy and an oxygen-sensitive nanoprobe.,” The Journal of neuroscience : the official journal of the Society for Neuroscience 35(8), 3663–3675 (2015).


## Instructions
1. Download VANmodel_0.0.zip from [here](https://drive.google.com/drive/u/1/folders/1ofeCST9HWSLrFGk8WpN9-0dzM-5mqY4R).
2. Unzip all the files to a folder and copy the folder name. For example, in the below case, the folder name is ‘C:\Users\xcheng17\Documents\VAN\’
![Figure](Example.png)
3.	Open the Matlab script VANmodelrun.m. Make sure you change ‘myfolder1’ to the folder you just saved your files. 
4.	If you are using Windows instead of Linux, make sure you change all the ‘/’ in all the scripts to ‘\’;
5.	You can specify the parameters you want, such as the magnetic field strength B0.
6.	Run the script. All the results, including intermediate results, are saved in the ‘data/mouse#/results/’ folder. Here the number ‘#’ is a number from 1 to 6 as you specified in the variable ‘mouseindex’.
7.	Do something else. The computation takes about 12-24 hours, depending on your computer speed.
8.	The gradient echo BOLD signal is the ‘BOLD_Gx’ variable in ‘MCBOLD_##.mat’.



