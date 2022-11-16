function [,] = Model_Calibration(P, Q_obs, ET_0, kc, K_sat, c, t_sub, t_sup, z, sw, s1, n, Q_b, A)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

cr = 1/1200 ; %cooling rate, table 3
n_iter = 100 ;



[Qsub, Qsup, R, I, s, L, ET] = hydro_model(P, Q_obs, ET_0, kc, K_sat, c, t_sub, t_sup, z, sw, s1, n, Q_b, A, phi,mode) ;

NS_old = 1-(sum((Q_obs-Qsub).^2)/sum((Q_obs-mean(Q_obs)).^2));

%sigma around 5\% of parameter range // sigma = walk coeff
sigma_K_sat = ;


for i = 1:n_iter

    %1. define a functional form of the temperature
    T_SA = t_sa(cr, i) ;

    %2. attribute arbitrary values to parameters
    
    K_sat_new = TruncNormRnd(K_sat, sigma_K_sat)

end



end