function [Qsub, Qsup, R, I, s, L, ET] = hydro_model(P, Q_obs, ET_0, kc, K_sat, c, t_sub, t_sup, z, sw, s1, n, Q_b, A, phi,mode)

%%% NE PAS MODIFIER CETTE FONCTION -> hydrological_model.m a la place :)

% This function computes the hydrological model (ie the discharge, runoff, infiltration,
%saturation, leaching and potential evapotranspiration. It can run the hydrological
%model for every hour (mode =0), every day in a month (mode =1), every month in a year (mode=2),
% and all years (mode = 3), depending on the mode given.
    
    Nyears=length(Q_obs)/(365*24);                   % number of years
    day_month=[31 28 31 30 31 30 31 31 30 31 30 31]; % number of days for each month
    hour_month=day_month*24;                         % number of hours for each month
    
    month_end=cumsum(day_month);                     % last day of each month
    month_start=month_end-day_month+1;               % first day of each month
    
    month_end_hour=cumsum(hour_month);               % last hour of each month   
    month_start_hour=month_end_hour-hour_month+1;    % fisrt hour of each month
    
    Q_hour_year=reshape(Q_obs,365*24,Nyears);
    P_hour_year=reshape(P,365*24,Nyears);

    if mode == 0   %model for every hour
       
        % Settings
        s = zeros(length(Q_obs),1);    % saturation vector
        s(1) = 0.5;                    % initial arbitrary value
        
        Vsup = zeros(length(Q_obs),1);  % superficial stored volume vector
        Vsup(1) = 0.5;                  % initial arbitrary value
        
        Vsub = zeros(length(Q_obs),1); %  sub-superficial stored volume vector
        Vsub(1) = 0.5;                  % initial arbitrary value
        
        N_hours_per_year = 365*24;  %number of hours per year
        ET_hour_year = zeros(N_hours_per_year,1); 
        L = zeros(length(Q_obs),1);
        
        % Infiltration
        K_sat_vect = ones(length(P),1)*K_sat;  % K_sat vector to compare with precipitation
        I = min(P,K_sat_vect);
        I_hour_year = reshape(I,N_hours_per_year,Nyears); %reshaping infiltration vector into matrix with hours in lines and years in columns
        
        % Runoff
        R = P-I;
        
        for y=1:Nyears
            for m=1:12
                for h=month_start_hour(m):month_end_hour(m)
                    t=(y-1)*N_hours_per_year+h; 
                    if s(t)<=sw && s(t)>=0
                        ET_hour_year(t)= 0;
                    elseif s(t)<=s1 && s(t)>sw
                        ET_hour_year(t)= kc(m)*ET_0(m)*((s(t)-sw)/(s1-sw));
                    elseif s(t)<=1 && s(t)>s1
                        ET_hour_year(t)= kc(m)*ET_0(m);
                    end
                    L(t)=K_sat*(s(t)^c);
                    s(t+1)=s(t)+((I(t)-ET_hour_year(t)-L(t))*(1/(n*z)));
                    Vsup(t+1)=Vsup(t)+R(t)-((t_sup^(-1))*Vsup(t));
                    Vsub(t+1)=Vsub(t)+L(t)-((t_sub^(-1))*Vsub(t));
                end
            end
        end
        
        
        ET=reshape(ET_hour_year,[],1);
        qsup=(t_sup^(-1)).*Vsup;
        qsub=(t_sub^(-1)).*Vsub;
        Qsup=A*qsup/3.6;
        Qsub=A*qsub/3.6;
    end
end
