function sig_k = sigma_k(a,sig_0,P,g_k,Gp)
    sig_k = (a*(sig_0^2))./(P*Gp*g_k);
end
