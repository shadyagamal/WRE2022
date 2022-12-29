function [Qout, V, l, Q_HU] = flood_control(Q_simulated, dt, Q_lim, Volume, level, Q_T, level_min, level_max, Cq_sluice, Cq_spill, L_s ,p)

% Constants
 
g = 9.81 ; % gravitational acceleration [m/s^2]


% Variables preallocation
V=zeros(length(Q_simulated),1);           % volume time series
l=zeros(length(Q_simulated),1);           % level time series
sluice_A=zeros(length(Q_simulated),1);    % sluice gate opening time series
Qout=zeros(length(Q_simulated),1);        % outflowing dicharge time series
Q_HU=zeros(length(Q_simulated),1);        % discharge for hydropower time series
Qg=zeros(length(Q_simulated),1);          % Sluice discharge


%%% Initial conditions %%%
V(1)=Volume(p+1);   %initial volume (level equal to the spillway level)
l(1)=p; %[m]
volume_max = Volume(level_max+1) %[m^3]


for t=1:dt:5000

    l(t+1) = level_volume(Volume,level,V(t));  %reservoir level


    %%  Computation of Q_HU [m^3/h] %%%

    ndays = fix(t/24);   %get the number of days elapsed
    

    pump_start = ndays*24 + 6;
    pump_stop = ndays*24 +22;
    
    if t>=pump_start && t<=pump_stop && l(t)>= level_min && l(t)<=level_max
        Q_HU(t) = Q_T; 
    else 
        Q_HU(t) = 0;
    end


    %%%  Sluice gate area computation %%%
    
    % computation of discharge through the gate 
    Q1_buffer = (V(t) + (Q_simulated(t)-Q_HU(t))*dt - volume_max)/dt;
    Q2_buffer = min(Q_lim,Q1_buffer);
    Qg(t)= max(Q347, Q2_buffer);  %[m^3/h]

    % computation of exit velocity applying energy equation between water
    % table in the reservoir and reservoir exit

    v_exit = sqrt(l(t)/(2*g))*3600; %exit velocity [m/h]
    
    % Sluice gate area [m^2]

    sluice_A(t) = Qg(t)/v_exit;


    %%% Total output discharge computation %%%

    if l(t) <= p 
        Qout(t) = Cq_sluice * sluice_A(t) * sqrt(2*g*l(t));
    elseif l(t)>p
        Qout(t) = Cq_sluice * sluice_A(t) * sqrt(2*g*l(t)) + Cq_spill * L_s * sqrt(2*g*(l(t)-p)^3);
    end
end
    

    %%% Storage equation integration %%%

    dV = (Q_simulated(t) - Qout(t) - Q_HU(t)) * dt;
    V(t+1) = V(t) + dV;