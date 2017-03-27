% Initialization for the project
% Add subdirectories to the searching path
% Authors: Yukun Chen, College of Information Science and Technology, Penn State University
% Contact: cykustc@gmail.com
os_str=computer;
% if strcmp(os_str,'MACI64')
gmmhmm_projectroot='.';
cd(gmmhmm_projectroot);
% elseif strcmp(os_str,'GLNXA64')
% cd('/gpfs/home/yzc147/work/GMMHMM');
% gmmhmm_projectroot='/gpfs/home/yzc147/work/GMMHMM';
% end
addpath(genpath('./utils'));
addpath(genpath('./exp'));
addpath(genpath('./src'));


%optim setting:
global optim_options 
optim_options = optimset('Display','off', 'LargeScale','off', 'Diagnostics','off');
global lpoptim_options 
lpoptim_options = optimset('Display','off', 'LargeScale','off', 'Diagnostics','off', 'Simplex', 'on', 'Algorithm','simplex');
global default_options 
default_options = optimset('Display','off', 'Diagnostics','off');
warning('off','optim:linprog:AlgOptsWillError')

addpath('/Users/yzc147/Documents/MATLAB/mosek/7/toolbox/r2012a/');

global problemMap problemSet
problemMap = containers.Map('KeyType', 'double', 'ValueType', 'double');
problemSet = {};