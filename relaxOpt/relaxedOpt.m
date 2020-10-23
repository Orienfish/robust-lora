% Call SNOPT for the relaxed version of the problem
setenv('SNOPT_LICENSE','~/Github/snopt-matlab/snopt7.lic');

% Important parameters
SF_cnt = 4;
CH_cnt = 8;
Ptx_cnt = 6;

M = 1;
sr_loc = csvread('sr_loc.csv');
gw_loc = csvread('gw_loc.csv');
sr_cnt = size(sr_loc, 1);
gw_cnt = size(gw_loc, 1);
c_ijk = zeros(sr_cnt, gw_cnt, SF_cnt);
for i = 0:SF_cnt-1
    f = 'cijk_' + string(i) + '.csv';
    d = csvread(f);
    c_ijk(1:end, 1:end, i+1) = d;
end

% Variable
x0 = zeros(gw_cnt, 1);

% Objective function
f = ones(1, gw_cnt);

% Linear inequality constraint: A * x <= b
A = - c_ijk(1:end, 1:end, SF_cnt);
b = - M * ones(sr_cnt, 1);

% Linear equality constraint: Aeq * x = beq
Aeq = [];
beq = [];

% Upper and lower bounds
lb = zeros(gw_cnt, 1);
ub = ones(gw_cnt, 1);

% Call the optimal solver
x = fmincon(@(x)(f*x), x0, A, b, Aeq, beq, lb, ub);
% [x,fval,INFO,output,lambda,states] = snsolve(@(x)(f*x), x0, A, b, Aeq, beq, lb, ub);