function d_hat = sense_target(s_t, s_q)
% S_q quads pos
% s_t target pos (the real one)
    d_s = norm(s_t-s_q);
    d_hat = d_s + sigma_k(d_s)*randn;
end

