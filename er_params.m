function params = er_params(nrows, ncols)
params.er_winsz = 7;
params.er_k = (params.er_winsz * 2 + 1)*2;
params.er_cnbins = 128;
params.tr_maxm = sqrt(nrows * nrows + ncols * ncols)/3;
params.er_sptol = 4;
params.er_cov_disturb = 1e-6;
params.er_cov_scale = 2;
params.er_ovlp_thre = 0.5;
params.er_diff_thre = 0.3;
params.er_sim_thre = 0.5;
params.er_szprior = 5e3;
params.er_score_thre = 0.2;
end