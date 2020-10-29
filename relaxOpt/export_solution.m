% Export the generated solution
function export_solution(x, sr_loc, gw_loc, params)
    % Export gateway placement
    gw_extract = [eye(params.gw_cnt), zeros(params.gw_cnt, params.var_cnt - params.gw_cnt)];
    gw_mask = logical(round(gw_extract * x));
    placed_gw_loc = gw_loc(gw_mask, :);
    fid = fopen('gw_relaxOpt.txt', 'w');
    for i = 1:size(placed_gw_loc, 1)
        fprintf(fid, '%f %f\n', placed_gw_loc(i, 1), placed_gw_loc(i, 2));
    end
    fclose(fid);

    % Export end device placement
    fid = fopen('sr_relaxOpt.txt', 'w');
    for i = 1:size(sr_loc, 1)
        % Find the assigned sf as the index of the max value
        sf_i = x(params.sf_st + (i-1) * params.SF_cnt + 1 : ...
            params.sf_st + i * params.SF_cnt);
        [val, sf] = max(sf_i);
        sf = sf - 1; % Change from starting from 1 to starting from 0
        % Find the assigned channel as the index of the max value
        %ch_i = x(params.ch_st + (i-1) * params.CH_cnt + 1 : ...
        %    params.ch_st + i * params.CH_cnt);
        %[val, ch] = max(ch_i);
        %ch = ch - 1; % Change from starting from 1 to starting from 0
        % Find the assigned tx power as the index of the max value
        tp_i = x(params.tp_st + (i-1) * params.Ptx_cnt + 1 : ...
            params.tp_st + i * params.Ptx_cnt);
        [val, tp] = max(tp_i);
        tp = params.Ptx_array(tp);
        
        fprintf(fid, '%f %f %d 20.0\n', sr_loc(i, 1), sr_loc(i, 2), ...
            sf);
    end
    fclose(fid);
end
