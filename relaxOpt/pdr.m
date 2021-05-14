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
hij = zeros(params.sr_cnt, params.gw_cnt);
for i = 1:params.sr_cnt
    % Get the current selections at i
    sf_i = x(params.sf_st + (i-1) * params.SF_cnt + 1 : ...
             params.sf_st + i * params.SF_cnt);
    ch_i = x(params.ch_st + (i-1) * params.CH_cnt + 1 : ...
             params.ch_st + i * params.CH_cnt);
    %tp_i = x(params.tp_st + (i-1) * params.TP_cnt + 1 : ...
    %         params.tp_st + i * params.TP_cnt);
    Nij = 0; % Number of nodes transmitting with the same SF and channel to gateway j
    for j = 1:params.gw_cnt
        for i_prime = 1:params.sr_cnt
            % Get the current CH selection at i prime
            ch_i_prime = x(params.ch_st + (i_prime-1) * params.CH_cnt + 1 : ...
                           params.ch_st + i_prime * params.CH_cnt);
            for k = 1:params.SF_cnt
                 % Get the current binary selection of SF k at i prime
                sf_i_prime_k = x(params.sf_st + (i_prime-1) * params.SF_cnt + k);
                for s = 1:params.TP_cnt
                    % Get the current binary selection of TP s at i prime
                    tp_i_prime_s = x(params.tp_st + (i_prime-1) * params.TP_cnt + s);
                    N_ij = N_ij + c_ijks(i_prime, j, k, s) 
                end
            end
        end
    end
    
    sf_i = x(params.sf_st + (i-1) * params.SF_cnt + 1 : ...
        params.sf_st + i * params.SF_cnt);
    
    
    hi(i) = 1 - exp(-2 * (params.Tk * sf_i) * Ni / params.Time);
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

c = params.PDRth - PDRi;
ceq = [];
end