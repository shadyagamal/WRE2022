function [,] = Model_Calibration(P, Q_obs, ET_0, kc, K_sat, c, t_sub, t_sup, z, sw, s1, n, Q_b, A)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here




[Qsub, Qsup, R, I, s, L, ET] = hydro_model(P, Q_obs, ET_0, kc, K_sat, c, t_sub, t_sup, z, sw, s1, n, Q_b, A, phi,mode) ;




end