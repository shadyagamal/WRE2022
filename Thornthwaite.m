function ET_0 = Thornthwaite(temperature, phi)
% Computes reference crop evapotranspiration ET0 
% (potential evapotranspiration) 
% through the Thornthwaite equation 
%   input is a time series of the temperatures


I= sum((temperature./5).^(1.514)) ; %heat index
a= 6.75*(10^(-7))*(I^3) - 7.71*(10^(-5))*(I^2) + 1.79*(10^(-2))*I + 0.49; %experimental exponent
D = 1:365;
delta=0.409*sin(2*pi*D/365-1.39);
w_s = acos(-tan(phi*pi/180)*tan(delta));
N_D = 24*w_s/pi; %number of daylight hours of day D
N_M =zeros(12,1); % mean daylight hours of month m

day_month=[31 28 31 30 31 30 31 31 30 31 30 31]; %"day_month": number of days for each month
month_end=cumsum(day_month);                     %"month_end": last day of each month
month_start=month_end-day_month+1;               %"month_start": first day of each month


for i=1:12
    N_M(i) = mean((N_D(month_start(i):month_end(i)))) ;
end

ET_0 = (16/12)*N_M .* ((10*temperature/I)).^a ./ day_month; % [mm/unit area]

end