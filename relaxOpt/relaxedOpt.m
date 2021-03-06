% Call SNOPT for the relaxed version of the problem
setenv('SNOPT_LICENSE','/home/xiaofan/Github/snopt-matlab/snopt7.lic');
addpath('~/Github/snopt-matlab/');
addpath('~/Github/snopt-matlab/util/');
clc;
clear;
close all;
warning('off','all');

% Important parameters
params.SF_cnt = 4;
params.CH_cnt = 8;
params.TP_cnt = 6;
params.Tk = [0.098, 0.175, 0.329, 0.616];
params.Time = 1200;
params.PDRth = 0.8;
params.M = 1;
params.Ptx_array = [5.0, 8.0, 11.0, 14.0, 17.0, 20.0];

sr_loc = csvread('sr_loc.csv');
gw_loc = csvread('gw_loc.csv');
params.sr_cnt = size(sr_loc, 1);
params.gw_cnt = size(gw_loc, 1);
c_ijks = zeros(params.sr_cnt, params.gw_cnt, params.SF_cnt, params.TP_cnt);
for i = 0:params.SF_cnt-1
    for j = 0:params.TP_cnt-1
        f = 'cijk_' + string(i) + '_' + string(j) + '.csv';
        d = csvread(f);
        c_ijks(1:end, 1:end, i+1, j+1) = d;
    end
end

% Variable
params.var_cnt = params.gw_cnt + params.sr_cnt * (params.SF_cnt + params.CH_cnt + params.TP_cnt);
gw_mask = [ones(1, params.gw_cnt, 1), zeros(1, params.var_cnt-params.gw_cnt)];
% (st, ed]
params.gw_st = 0; params.gw_ed = params.gw_cnt;
params.sf_st = params.gw_ed; params.sf_ed = params.sf_st + params.SF_cnt * params.sr_cnt;
params.ch_st = params.sf_ed; params.ch_ed = params.ch_st + params.CH_cnt * params.sr_cnt;
params.tp_st = params.ch_ed; params.tp_ed = params.tp_st + params.TP_cnt * params.sr_cnt;
x0 = zeros(params.var_cnt, 1);

% Objective function
f = gw_mask;

% Linear inequality constraint: A * x <= b
A1 = - [c_ijks(1:end, 1:end, params.SF_cnt, params.TP_cnt), ...
    zeros(params.sr_cnt, params.var_cnt-params.gw_cnt)];
b1 = - params.M * ones(params.sr_cnt, 1);
[A2, b2] = lifetimeConstraint(params);
A = [A1; A2];
b = [b1; b2];

% Linear equality constraint: Aeq * x = beq
[Aeq, beq] = validConstraint(params);

% Upper and lower bounds
lb = zeros(params.var_cnt, 1);
ub = ones(params.var_cnt, 1);

% Call the optimal solver
tic
%x = fmincon(@(x)(f*x), x0, A, b, Aeq, beq, lb, ub, @(x)pdr(x, c_ijk, params));
[x,fval,INFO,output,lambda,states] = snsolve(@(x)(f*x), x0, A, b, ...
    Aeq, beq, lb, ub, @(x)pdr(x, c_ijks, params));
exeTime = toc;
%x

% Print the results
gw_extract = [eye(params.gw_cnt), zeros(params.gw_cnt, params.var_cnt - params.gw_cnt)];
gw_mask = logical(round(gw_extract * x));
fid = fopen('result.txt', 'a+');
fprintf(fid, '%f,%d,%f\n', f*x, sum(gw_mask), exeTime);
fclose(fid);
export_solution(x, sr_loc, gw_loc, params);
plot_solution(sr_loc, gw_loc(gw_mask, 1:end));

% Plot the solution in the grid space
function plot_solution(sr_loc, gw_loc)
    % intialization
    sr_cnt = size(sr_loc, 1);  % number of end devices
    gw_cnt = size(gw_loc, 1);  % number of gateways

    % scatter
    color_map = [[0 0.4470 0.7410]; ...       % blue
        [0.8500 0.3250 0.0980]; ...           % orange
        [0.3010 0.7450 0.9330]; ...           % ryan
        [0.4660 0.6740 0.1880]; ...           % green
        [0.6350 0.0780 0.1840]];              % red
    nodes = vertcat(sr_loc, gw_loc);           % append PoIs
    sz_nodes = repmat(80, size(nodes, 1), 1); % const size for nodes
    %color_idx = vertcat(sol.x+sol.s, repmat(3, size(O, 1), 1), 4) + 1;
    color_idx = vertcat(ones(sr_cnt, 1), repmat(2, gw_cnt, 1));
    color_idx = round(color_idx);
    color_nodes = vertcat(color_map(color_idx, :));
    figure;
    scatter(nodes(:, 1), nodes(:, 2), sz_nodes, color_nodes, 'filled', ...
        'LineWidth', 2);
end

% Get the matrix for lifetime constraint
function [A2, b2] = lifetimeConstraint(params)
A2 = zeros(params.sr_cnt, params.var_cnt);
for i = 1:params.sr_cnt
    % Converted lifetime constraint: sf_i^4+tp_i^5+tp_i^6 <= 1
    A2(i, params.sf_st+(i-1)*params.SF_cnt+4) = 1; % sf_i^4
    A2(i, params.tp_st+(i-1)*params.TP_cnt+5) = 1; % tp_i^5
    A2(i, params.tp_st+(i-1)*params.TP_cnt+6) = 1; % tp_i^6
end
b2 = ones(size(A2, 1), 1);
end

% Get the matrix for validation constraint
function [Aeq, beq] = validConstraint(params)
Aeq = zeros(3*params.sr_cnt, params.var_cnt);
for i = 1:params.sr_cnt
    % SF
    for k = 1:params.SF_cnt
        Aeq((i-1)*3+1, params.sf_st+(i-1)*params.SF_cnt+k) = 1; % sf_i^{1-4}
    end
    % CH
    for q = 1:params.CH_cnt
        Aeq((i-1)*3+2, params.ch_st+(i-1)*params.CH_cnt+q) = 1; % ch_i^{1-8}
    end
    % TP
    for s = 1:params.TP_cnt
        Aeq((i-1)*3+3, params.tp_st+(i-1)*params.TP_cnt+s) = 1; % ch_i^{1-8}
    end
end
beq = ones(size(Aeq, 1), 1);
end