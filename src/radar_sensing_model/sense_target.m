function d_hat = sense_target(s_t, S_q)
% S_q quads pos
% s_t target pos (the real one)
    d_s = norms(s_t-S_q,2);
    d_hat = d_s + sigma_k(d_s).*randn(1,length(d_s));
end

