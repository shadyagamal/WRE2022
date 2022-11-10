function T_SA = t_sa(cr,i)
%t_sa Functional form for the temperature
%   Where i counts the iterations of the calibration procedure, while cr is a cooling rate
T_SA = exp(-cr*i);
end

