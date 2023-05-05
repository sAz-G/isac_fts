function [idx, idy,min_val] = get_min(D_est,x_t,y_t,x_jt,y_jt,params)
    [idx, idy,min_val] = get_min_gridsearch(D_est,x_t,y_t,x_jt,y_jt,params);
end

