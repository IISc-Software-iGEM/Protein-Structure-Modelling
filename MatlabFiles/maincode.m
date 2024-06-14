% This program calls buildprot
%
% Babak Alipanahi
% University of Waterloo
% May 12, 2010
% eidt: -June 16, 2010

clear all;
clc;
close all;


%addpath('helperfunctions');
%addpath('inputreaders');
%addpath('sdpstuff');
%addpath('refinement');
%addpath('lib');
%addpath('../packages/hanso');

cd ..
addpath(genpath(pwd));
cd code;

load '../data/ainfo.mat'

% Drop side-chain hydrogen atoms or not
%hydrogen_omission = 1;
hydrogen_omission = 0;

% Number of BFGS iterations
gd_tol  = 10^-9;
cg_iter = 500;
f = [10 10 10 10 10]; % [f_hb, f_tau, f_tal, f_vdw, f_tas]



fprintf('==========================================================================\n')
fprintf('------------SPROS: SDP-based Protein Structure Determination--------------\n')



tstart = tic;

%protein_name  = 'ubiq';
%protein_name  = 'talin';
%protein_name = 'murine';
%protein_name = 'np';


in_hbond_file = '';
in_max_res = [];
in_min_res = [];

protein_name = '1qhk';
in_seq_file = '1qhk_23_29.seq'; 
in_upl_file = {'1qhk_concat_dist_23_29.upl'}; 
in_ang_file = ''; 
in_hbond_file = '1qhk_concat_hBond_23_29.upl'; 

%in_seq_file = '1qhk.seq'; 
%in_upl_file = {'1qhk_concat_dist.upl'}; 
%in_ang_file = ''; 
%in_hbond_file = '1qhk_concat_hBond.upl'; 

%in_seq_file = '1qhk_21_31.seq'; 
%in_upl_file = {'1qhk_concat_dist_21_31.upl'}; 
%in_ang_file = ''; 
%in_hbond_file = '1qhk_concat_hBond_21_31.upl'; 

protein_path = ['../proteins/' protein_name '/'];
fprintf('*************************************************************************\n');
fprintf('Protein: %s\n', protein_name);
fprintf('-Reading input files...\n')
% Reading input data
%==========================================================================
seq_file = [protein_path in_seq_file];
[seq, num] = seq_reader(seq_file);

%print seq and num
seq
num