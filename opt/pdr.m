function [c, ceq] = pdr(x, c_ijks, params)
% Evaluate PDR
% 
% Args:
%   x: variables including gateway placement and SF, channel and tx power
%      configuration
%   c_ijks: a binary matrix showing feasibility of i reaching j with SF k
%           and Tx Power s
%   params: important parameters
%
% Return:
%   c, ceq: nonlinear constaints, inequality and equality

% Split various variables from x
gw_extract = [eye(params.gw_cnt), zeros(params.gw_cnt, params.var_cnt - params.gw_cnt)];
g = gw_extract * x;

% Compute collision probability at end device i 
hi = zeros(params.sr_cnt, 1);
for i = 1:params.sr_cnt
    Ni = 0; % Number of nodes transmitting with the same SF and channel
    for j = 1:params.sr_cnt
        if j == i % j != i in calculating collision probabilities
            continue;
        end
        
        % Traverse all possible combinations of SF and channel
        for k = 1:params.SF_cnt
            for q = 1:params.CH_cnt
                sf_ik = x(params.sf_st + (i-1) * params.SF_cnt + k);
                sf_jk = x(params.sf_st + (j-1) * params.SF_cnt + k);
                ch_iq = x(params.ch_st + (i-1) * params.CH_cnt + q);
                ch_jq = x(params.ch_st + (j-1) * params.CH_cnt + q);
                Ni = Ni + sf_ik * sf_jk * ch_iq * ch_jq;
            end
        end
    end
    % Get the current SF selection at i
    sf_i = x(params.sf_st + (i-1) * params.SF_cnt + 1 : ...
        params.sf_st + i * params.SF_cnt);
    hi(i) = 1 - exp(-2 * (params.T_k * sf_i) * Ni / params.Time);
    %for k = 1:params.SF_cnt
    %    fprintf('%f ', sf_i(k));
    %end
    %fprintf('Collision probability: %f\n', hi);
end

% Compute transmission reliability from i to j
% combining collision probability hi and reachability c_ijks
Pij = zeros(params.sr_cnt, params.gw_cnt);
for i = 1:params.sr_cnt
    sf_i = x(params.sf_st + (i-1) * params.SF_cnt + 1 : ...
        params.sf_st + i * params.SF_cnt);
    tp_i = x(params.tp_st + (i-1) * params.TP_cnt + 1 : ...
        params.tp_st + i * params.TP_cnt);
    % Get a SF_cnt * TP_cnt matrix showing the result of every pair
    sfk_tps = sf_i * reshape(tp_i, [1, params.TP_cnt]);
    for j = 1:params.gw_cnt
        c_ij = c_ijks(i, j, 1:end, 1:end); % c_ij is SF_cnt * TP_cnt
        c_ij = squeeze(c_ij); % Remove dimension of length 1
        Pij(i, j) = sum(sfk_tps .* c_ij, 'all') * (1 - hi(i));
    end
end

% Compute PDR at i
PDRi = zeros(params.sr_cnt, 1);
for i = 1:params.sr_cnt
    P = 1 - Pij(i, 1:end)' .* g;
    PDRi(i) = 1 - prod(P);
end
PDRi'

c = params.PDR_th - PDRi;
ceq = [];
end