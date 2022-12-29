function P = P_prod(Q_pipe, D_pipe, L_pipe, ks, delta_z_empty, l, eta, )
 
% Constants
nu = 1.79 * 10^-5;   % water kinematic viscosity [m/s^2]
rho = 1000;   %water density [kg/m^3] 
g = 9.81  %gravitation acceleration [m/s^2]


V_pipe = (Q_pipe/3600)/(pi*(D_pipe/2)^2); % [m/s] water velocity in the pipe

Re = V_pipe * D_pipe / nu;  % Reynold's number [-]

if Re<= 2000  % laminar flow
   
    f= Re/64;

elseif Re>2000  % turbulent flow

    f = 0.25/(log10((ks/(3.7*D_pipe))+5.74/Re^0.9)^2)

end

delta_hl = f * (L_pipe/D_pipe)*(V_pipe^2/(2*g)); % linear head loss (Darcy-Weisbach)

 
delta_hs = 1.5 * V_pipe^2 / (2*g);  % singular head loss

delta_h = delta_hs + delta_hl;  % total head loss

h_tot = delta_z_empty + l;  %[m] initial total head

h_t = h_tot - delta_h;  % net head of the turbine [m]

 
P_turbine_J = eta * rho * g .* Q_HU /3600 .* h_t ;  % [J] hydropower production

nyears_P = length(P_turbine_J)/(365*24);  %number of years for hydropower production

P_years = zeros(nyears_P,1);

for y = 0:nyears_P-1
    P_years(y+1) = sum(P_turbine_J(y*24*365+1:(y+1)*24*365));  % [J/year]
    P_years(y+1) = P_years(y+1)/(365*24*10^9);    %[GWh]
end

P = mean(P_years);
