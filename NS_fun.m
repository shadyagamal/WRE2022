
function NS = NS_fun(Q_obs,Q_mod)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
NS = 1- ((sum((Q_obs-Q_mod).^2))/(sum((Q_obs-mean(Q_obs)).^2)));
end
