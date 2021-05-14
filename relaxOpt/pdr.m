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

% Compute collision probability transmitting from end device i to gateway j
h_ij = zeros(params.sr_cnt, params.gw_cnt);
for i = 1:params.sr_cnt
    % Get the current selections at i
    sf_i = x(params.sf_st + (i-1) * params.SF_cnt + 1 : ...
             params.sf_st + i * params.SF_cnt);
    ch_i = x(params.ch_st + (i-1) * params.CH_cnt + 1 : ...
             params.ch_st + i * params.CH_cnt);
    %tp_i = x(params.tp_st + (i-1) * params.TP_cnt + 1 : ...
    %         params.tp_st + i * params.TP_cnt);
    for j = 1:params.gw_cnt
        % Init the number of nodes transmitting with the same SF and channel 
        % to gateway j, considering feasibility by integrating c_ijks
        N_ij = 0;
        for i_prime = 1:params.sr_cnt
            % Get the current selections at i prime
            sf_i_prime = x(params.sf_st + (i_prime-1) * params.SF_cnt + 1 : ...
                           params.sf_st + i_prime * params.SF_cnt);
            ch_i_prime = x(params.ch_st + (i_prime-1) * params.CH_cnt + 1 : ...
                           params.ch_st + i_prime * params.CH_cnt);
            tp_i_prime = x(params.tp_st + (i_prime-1) * params.TP_cnt + 1 : ...
                           params.tp_st + i_prime * params.TP_cnt);
            % Get one number for sf_i = sf_i_prime
            SF_equal = sf_i * sf_i_prime';
            % Get one number for ch_i = ch_i_prime
            CH_equal = ch_i * ch_i_prime';
            % Accumulate N_ij
            N_ij = N_ij + SF_equal * c_ijks(i_prime, j, :, :) * ...
                tp_i_prime' * CH_equal; 
        end
        h_ij(i, j) = exp(-2 * (params.Tk * sf_i) * N_ij / params.Time);
    end
end

fprintf('Non-Collision Probability: %f\n', h_ij);

% Compute transmission reliability from i to j
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

c = params.PDRth - PDRi;
ceq = [];
end