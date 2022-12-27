function [Q_mod,R,I, s, ET, L] = hydrological_model(P, ET_0, kc, K_sat, c, t_sub, t_sup, z, sw, s1, n, Q_b, A,doTest)

% This function computes the hydrological model (ie the discharge, runoff, infiltration,
%saturation, leaching and potential evapotranspiration in [m^3/s]). It can run the hydrological
%model for every hour (mode =0), every day in a month (mode =1), every month in a year (mode=2),
% and all years (mode = 3), depending on the mode given.
% ALL UNITS MUST BE IN m and h
  
    Nyears=length(P)/(365*24);                       % number of years
    N_hours_per_year = 365*24;                       % number of hours per year 
    day_month=[31 28 31 30 31 30 31 31 30 31 30 31]; % number of days for each month
    hour_month=day_month*24;                         % number of hours for each month

    %month_end=cumsum(day_month);                     % last day of each month
    %month_start=month_end-day_month+1;               % first day of each month

    month_end_hour=cumsum(hour_month);               % last hour of each month   
    month_start_hour=month_end_hour-hour_month+1;    % fisrt hour of each month



    % Settings
    s = zeros(size(P));    % saturation vector
    s(1) = 0.5;            % initial arbitrary value

    Q_mod = zeros(size(P));    %Discharge
    R = zeros(size(P));    %Run off
    I = zeros(size(P));    %Infiltration
    L = zeros(size(P));    %Leaching
    ET = zeros(size(P));   %Evapotranspiration
    qsup = zeros(size(P)); %Superficial specific discharge
    qsub = zeros(size(P)); %Sub-superficial specific discharge
    Vsup = zeros(size(P)); %Superficial storage
    Vsub = zeros(size(P)); %Sub-superficial storage
    Qsup = zeros(size(P)); %Superficial discharge
    Qsub = zeros(size(P)); %Sub-superficial discharge
    dt = 1;



    for y=1:Nyears
        for m=1:12
            for h=month_start_hour(m):month_end_hour(m)
                t=(y-1)*N_hours_per_year+h;

                I(t)=min(P(t),K_sat); %m/h Infiltration
                R(t)=P(t)-I(t); %m/h Runoff

                if s(t)<=sw && s(t)>=0
                    ET(t)= 0; %m/h
                elseif s(t)<=s1 && s(t)>sw
                    ET(t)= kc(m)*ET_0(m)*((s(t)-sw)/(s1-sw)); %m/h
                elseif s(t)<=1 && s(t)>s1
                    ET(t)= kc(m)*ET_0(m); %m/h
                elseif s(t)>1
                        disp('Warning : s > 1')
                elseif s(t)<0
                        disp('Warning : s < 0')
                end
                L(t)=K_sat*(s(t)^c); %m/h
                s(t+dt)=s(t)+(I(t)-ET(t)-L(t))*dt/(n*z);

                qsup(t)=Vsup(t)/t_sup;
                qsub(t)=Vsub(t)/t_sub;

                Vsup(t+dt)=Vsup(t)+(R(t)-qsup(t))*dt; %m3
                Vsub(t+dt)=Vsub(t)+(L(t)-qsub(t))*dt; %m3

                Qsup(t)=A*qsup(t); % m3/h
                Qsub(t)=A*qsub(t); % m3/h
                Q_mod(t) = Qsub(t)+Qsup(t)+Q_b; % m3/h
            end
        end
    end

 %compute total discharge Q    
if doTest     % if doTest=1, check mass balance 

    P_tot=sum(P)*dt;
    R_tot=sum(R)*dt;
    L_tot=sum(L)*dt;
    ET_tot=sum(ET)*dt;

    %"testS" balance for the root zone (input/output). testS close to unity: necessary (but not sufficient) 
    %condition for the implementation to be correct 
    testS=P_tot/(ET_tot+R_tot+L_tot+n*z*(s(end)-s(1))) ;

    %"testQ" balance for the whole system (input/output). testQ close to unity: necessary (but not sufficient) 
    %condition for the implementation to be correct 
    testQ=sum(P-ET)/(sum(qsup+qsub)+n*z*(s(end)-s(1))+Vsup(end)+Vsub(end)) ;
end   

end