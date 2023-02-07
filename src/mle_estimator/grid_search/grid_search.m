% a script to test the functionality of the mle estimator. functionality still has to be
% tested.

clear, clc;
D_est = [1,2,3,4,6,7,8,4].*randn(1,8);
[x_t,y_t] = grid_search(1500,1500,100,100);
x_jt= [1,2,3,4,0,5,7,8].*randn(1,8);
y_jt= [1,2,3,4,6,5,4,6].*randn(1,8);
H = 100;
[x,y] = get_min(D_est,x_t,y_t,x_jt,y_jt, H);