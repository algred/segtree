function params = er_params(nrows, ncols)

params.er_winsz = 7;
params.er_k = (params.er_winsz * 2 + 1)*2;
params.er_cnbins = 128;
params.tr_maxm = sqrt(nrows * nrows + ncols * ncols)/3;
params.er_sptol = 4;
params.er_cov_disturb = 1e-6;
params.er_cov_scale = 1;
params.er_ovlp_thre = 0.5;
params.er_diff_thre = 0.4;
params.er_sim_thre = 0.5;
params.er_szprior = 5e3;
params.er_score_thre = 0.2;
params.er_delta = 3;
params.er_match_thre = 0.7;
params.er_min_len_ratio = 0.5;
%params.er_sz_chg_thre = 1.5;
params.minsz_ratio = 0.003;
params.min_motion = 0.8;
params.maxsz_ratio = 0.25;
params.min_shrink_ratio = 0.67;
params.DISPLAY = 0;
params.fg_ratio = 0.3;
params.bw_quant_num = 5;
params.bw_sigma_blur = 30;
params.pr_maxbgm = 30;
params.pr_cnbins = 128;
params.verbose = 1;

end