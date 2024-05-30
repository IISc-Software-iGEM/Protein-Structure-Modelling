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
max_res = max(num);
min_res = min(num);
    
num_upl = size(in_upl_file, 2);
upl_file = cell(1, num_upl);
for i = 1:num_upl
    upl_file{i} = [protein_path in_upl_file{i}];
end
if ~isempty(in_hbond_file)
    hbond_file = [protein_path in_hbond_file];
    hbond_write_file = [protein_path protein_name '_hbo.upl'];
    %hbond_reader(hbond_file,hbond_write_file);
    hbond_reader_me(hbond_file,hbond_write_file);
    num_upl = num_upl + 1;
    upl_file{num_upl} = hbond_write_file;
end
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
raw_up    = cell(1, num_upl);
raw_up_ho = cell(1, num_upl);
temp_min_res = +inf;
temp_max_res = -inf;
for i = 1:num_upl
    raw_up{i} = dist_reader(upl_file{i}, A);
    if hydrogen_omission
        raw_up_ho{i} = dist_reader(upl_file{i}, A, hydrogen_omission);
    end
    temp_min_res = min(temp_min_res, min([raw_up{i}.tres raw_up{i}.sres]));
    temp_max_res = max(temp_max_res, max([raw_up{i}.tres raw_up{i}.sres]));
end
if isempty(in_min_res)
    min_res = max(temp_min_res-1, min_res);
else
    min_res = in_min_res;
end
if isempty(in_max_res)
    max_res = min(temp_max_res+1, max_res);
else
    max_res = in_max_res;
end
    
% Remove informationless (w/o any constraints)
% parts from and N- and C-terminus
ind_del_N = num < min_res;
ind_del_C = num > max_res;
ind_del   = ind_del_N | ind_del_C;
seq(ind_del) = [];
num(ind_del) = [];
    
% dihedral angle constraints
ang_file = [protein_path in_ang_file];
[phi_cons, psi_cons] = ang_reader(ang_file, num);
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
pdb_write_file = [protein_path protein_name '.pdb'];
%==========================================================================
tread = toc(tstart);
fprintf('\tdone: %4.1f sec\n', tread)

fprintf('-Sampling a random molecule...\n')
% Generating a random structure
%==========================================================================
[phi, psi] = ang_sampler(seq,phi_cons,psi_cons);
switch hydrogen_omission
    case 0
        [rand_X, Comp] = ibuildprot(seq,num,phi,psi,A);
        wh_rand_X = rand_X;
        wh_Comp   = Comp;
    case 1
        [wh_rand_X, wh_Comp, rand_X, Comp] = ibuildprot(seq,num,phi,psi,A);
end

if exist(ang_file, 'file')
    [ang_lo_cons, ang_up_cons] = ang_dist_conmaker(phi_cons,psi_cons,wh_Comp);
end
%==========================================================================
trand = toc(tstart);
fprintf('\tdone: %4.1f sec\n',trand - tread)

fprintf('-Reducing the SDP problem...\n')

% Perform the reduction
%==========================================================================
switch hydrogen_omission
    case 0
        [U, Comp.cliq_dims] = reducer(rand_X,Comp);
        wh_Comp = Comp;
    case 1
        dont_compute_U = 1;
        [U, Comp.cliq_dims] = reducer(rand_X,Comp);
        [wh_U, wh_Comp.cliq_dims] = reducer(wh_rand_X,wh_Comp,dont_compute_U);
end
%==========================================================================
treduc = toc(tstart);
fprintf('\tdone: %4.1f sec\n',treduc - trand)
fprintf('-Forming constraints...\n')
% Generating upper and lower bounds constraints
%==========================================================================
start = 0;
wh_up_bounds = nan(50000,4);
for i = 1:num_upl
    temp_upl = upper_maker(raw_up{i}, wh_Comp);
    wh_up_bounds(start+1:start+size(temp_upl,1),:) = temp_upl;
    start = start + size(temp_upl,1);
end
wh_up_bounds(isnan(wh_up_bounds(:,1)), :) = [];

if hydrogen_omission
    start = 0;
    ho_up_bounds = nan(50000,4);
    for i = 1:num_upl
        temp_upl = upper_maker(raw_up_ho{i}, wh_Comp);
        ho_up_bounds(start+1:start+size(temp_upl,1),:) = temp_upl;
        start = start + size(temp_upl,1);
    end
    ho_up_bounds(isnan(ho_up_bounds(:,1)), :) = [];
end

% adding torsion-angle constraints
if exist(ang_file, 'file')
    wh_up_bounds = [wh_up_bounds; ang_up_cons];
    wh_sdp_lo_bounds = ang_lo_cons;
end
    
switch hydrogen_omission
    case 0
        % equality cons
        wh_eq_cons = equality_con_former(rand_X,Comp);
        eq_cons = wh_eq_cons;
        % upper bounds
        up_bounds = wh_up_bounds;
        % vdw bounds
        wh_vdw_bounds = vdw_bound_maker(wh_Comp);
        vdw_bounds    = wh_vdw_bounds;
        sdp_lo_bounds = wh_sdp_lo_bounds;
    case 1
        % equality cons
        wh_eq_cons = equality_con_former(wh_rand_X,wh_Comp);
        eq_cons    = equality_con_former(rand_X,Comp);
        % upper bounds
        up_bounds = map_bounds(ho_up_bounds,Comp.atoms_map);
        if exist(ang_file,'file')
            ang_up_bounds = map_bounds(ang_up_cons,Comp.atoms_map);
            up_bounds = [up_bounds; ang_up_bounds];
        end
        % vdw bounds
        wh_vdw_bounds = vdw_bound_maker(wh_Comp);
        vdw_bounds    = vdw_bound_maker(Comp);
        sdp_lo_bounds = map_bounds(ang_lo_cons,Comp.atoms_map);
        % a bug that need to be fixed sometimes
        % when TYR is HH it is mapped to its OH
        up_bounds(:,4) = wh_up_bounds(:,4);
end
    
lo_bounds     = vdw_bounds;
wh_lo_bounds  = wh_vdw_bounds;
if ~exist('sdp_lo_bounds', 'var')
    sdp_lo_bounds = [];
else
    lo_bounds    = [lo_bounds; sdp_lo_bounds];
    wh_lo_bounds = [wh_lo_bounds; wh_sdp_lo_bounds];
end
%==========================================================================
tbounds = toc(tstart);
fprintf('\tdone: %4.1f sec\n',tbounds - treduc)


fprintf('-Solving the SDP problem...\n')
fprintf('\n');
% Solve the SDP
%==========================================================================
%==========================================================================
[rawX, eV, eS] = solve_sdpt3(U,eq_cons,up_bounds,sdp_lo_bounds,[],f);
%--[rawX_cvx_spros, eigens_cvx_spros, slacks_cvx_spros] = solve_by_cvx_SPROS(U, eq_cons, up_bounds, sdp_lo_bounds, [], f);
%--[rawX_cvx, eigens_cvx, slacks_cvx] = solve_by_cvx_rank(size(U,1), eq_cons, up_bounds, sdp_lo_bounds, [], f);
%==========================================================================
%figure(1)
%bspros = bar(eV);

%figure(2)
%b_cvx_rank_spros = bar(eigens_cvx_spros);

%figure(3)
%b_cvx = bar(eigens_cvx);
%==========================================================================
%orig_rawX = rawX_cvx_spros;
%rawX = rawX_cvx_spros(1:3,:);

orig_rawX = rawX;
rawX = rawX(1:3,:);
%rawX = rawX_cvx(1:3,:);
%rawX = rawX_cvx_spros(1:3,:);
%==========================================================================
tsdp = toc(tstart);
fprintf('\tdone: %4.1f sec\n',tsdp - tbounds)
fprintf('\n\n')
disp('*********************************************************************')
% Analysis of the output
%==========================================================================
fprintf('Violations (raw)\n')
check_eq_cons = equality_con_former(rand_X, Comp, 2);
report = protchecker(rawX(1:3,:),Comp,check_eq_cons,lo_bounds,up_bounds,1);
disp('*********************************************************************')
%==========================================================================
% Post-processing
%==========================================================================
fprintf('Post-Processing...\n')
fprintf('\n');
switch hydrogen_omission
    case 0
        W = [2 1 1 -1];
    case 1
        W = [2 1 1 -1];
end
% Phase I - GD
% Refinement by HANSO
[pX, hinfo] = hanso_post_processing(rand_X, rawX, Comp, lo_bounds, up_bounds, W, f);
fprintf('Violations (GD-I)\n')
p_report = protchecker(pX,Comp,check_eq_cons,lo_bounds,up_bounds,1);
%---------------------------------------------------------------------------
% simple fix using the fact that most residues lie on the left half of
% Ramachandran plot
if sum(p_report.phi(~isnan(p_report.phi)) > 0) > 0.5*length(p_report.phi)
    pX(1,:) = -pX(1,:);
end
p_report = protchecker(pX,Comp,check_eq_cons,lo_bounds,up_bounds,0);
fprintf('\nCorrecting chiralities:\n')
pXc = chirality_correction(pX,Comp,p_report.chiral);
fprintf('Violations (after fixing chiralities)\n')
p_report = protchecker(pXc,Comp,check_eq_cons,lo_bounds,up_bounds,1);
disp('*********************************************************************')

[pXc, c_info] = hanso_post_processing(rand_X,pXc,Comp,lo_bounds,up_bounds,W,f);
fprintf('Violations (GD-II)\n')
pc_report = protchecker(pXc,Comp,check_eq_cons,lo_bounds,up_bounds,1);
fprintf('\nCorrecting chiralities:\n')
pXc = chirality_correction(pXc,Comp,pc_report.chiral);
disp('*********************************************************************')

    
if hydrogen_omission
    disp('*********************************************************************')
    fprintf('Putting back hydrogen atoms...\n')
    pX_wh = hydrogen_mapper(pXc,wh_rand_X,Comp,wh_Comp,A);
    check_eq_cons_wh = equality_con_former(wh_rand_X,wh_Comp,2);
    fprintf('Violations (After putting back hydrogens)\n')
    f_report = protchecker(pX_wh,wh_Comp,check_eq_cons_wh,wh_lo_bounds,wh_up_bounds,1);
    fprintf('\n');
    [fX, wh_info] = hanso_post_processing(wh_rand_X,pX_wh,wh_Comp,wh_lo_bounds,wh_up_bounds,W,f);
    fprintf('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
    fprintf('\nViolations (FINAL)\n')
    fprintf('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
    final_report = protchecker(fX,wh_Comp,check_eq_cons_wh,wh_lo_bounds,wh_up_bounds,1);
    disp('*********************************************************************')
end
%==========================================================================
sc_flag = 1;
switch hydrogen_omission
    case 0
        pdb_writer(pdb_write_file,pXc,Comp)
    case 1
        pdb_writer(pdb_write_file,fX,wh_Comp)
end
fprintf('\tOverall time: %4.1f sec\n',toc(tstart))
fprintf('done!\n')
fprintf('==========================================================================\n')
    


