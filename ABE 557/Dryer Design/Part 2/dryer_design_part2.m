%Clear command window, clear variable workspace, and close all figures
clc
clear all
close all

%Semolina Constants%
M0s = [11.8, 6.45, 5.92, 3.53]; % Data on monolayer constant
Ks = [0.65, 0.71, 0.72, 0.76]; % Data on K constant
Cs = [4.21, 8.77, 10.1, 200.1]; % Data on C constant


%Farina Constants
M0f = [9.16, 6.01, 5.23, 4.76];
Kf = [0.65, 0.69, 0.71, 0.73];
Cf = [14.06, 31.95, 37.04, 137.0];

%Create temperature Vector
T_C = [20, 35, 50, 60]; % Temps in C
T_K = T_C + 273; % Temps in Kelvin

%Define generic constants
R = 8.314; %J/mol-K

%Apply a linear regression analysis to each data set against the inverse of
%the temperature
coeffs_s_M0 = polyfit(1./T_K,log(M0s),1);
coeffs_s_K = polyfit(1./T_K,log(Ks),1);
coeffs_s_C = polyfit(1./T_K,log(Cs),1);
coeffs_f_M0 = polyfit(1./T_K,log(M0f),1);
coeffs_f_K = polyfit(1./T_K,log(Kf),1);
coeffs_f_C = polyfit(1./T_K,log(Cf),1);

%Calculate Parameters for semolina%
%See the backround section of the report
M0o_s = exp(coeffs_s_M0(2)); 
Ko_s = exp(coeffs_s_K(2));
Co_s = exp(coeffs_s_C(2));
H_M0_s = R*coeffs_s_M0(1);
H_K_s = R*coeffs_s_K(1);
H_C_s = R*coeffs_s_C(1);

%Calculate Patameters for ferina
M0o_f = exp(coeffs_f_M0(2));
Ko_f = exp(coeffs_f_K(2));
Co_f = exp(coeffs_f_C(2));
H_M0_f = R*coeffs_f_M0(1);
H_K_f = R*coeffs_f_K(1);
H_C_f = R*coeffs_f_C(1);

%Define the temperature vector to plot against
T0 = 0;
Tf = 100;
temp_profile = zeros(1,Tf); 

%populate the vector with the necessary temperatures in celsius
for i = 1:1:length(temp_profile)-1
    temp_profile(i+1) = temp_profile(i)+1; %Celcius
end

%Initialize the vectors for the GAB model constants as a function
%of temperature in our system.
M0_s_vect = zeros(1,length(temp_profile));
K_s_vect = zeros(1,length(temp_profile));
C_s_vect = zeros(1,length(temp_profile));
M0_f_vect = zeros(1,length(temp_profile));
K_f_vect = zeros(1,length(temp_profile));
C_f_vect = zeros(1,length(temp_profile));

%Iterate through the temperatures and calculate the GAB model equation
%constnats that are associated with each temp.

%Semoline m0
for i = 1:1:length(M0_s_vect);
    T = temp_profile(i);
    M0_s_vect(i) = para_calc(M0o_s,H_M0_s,T+273);
end

%Semolina K
for i = 1:1:length(K_s_vect);
    T = temp_profile(i);
    K_s_vect(i) = para_calc(Ko_s,H_K_s,T+273);
end

%Semolina C
for i = 1:1:length(C_s_vect);
    T = temp_profile(i);
    C_s_vect(i) = para_calc(Co_s,H_C_s,T+273);
end

%Farina M0
for i = 1:1:length(M0_f_vect);
    T = temp_profile(i);
    M0_f_vect(i) = para_calc(M0o_f,H_M0_f,T+273);
end

%Farina K
for i = 1:1:length(K_f_vect);
    T = temp_profile(i);
    K_f_vect(i) = para_calc(Ko_f,H_K_f,T+273);
end

%Farina C
for i = 1:1:length(C_f_vect);
    T = temp_profile(i);
    C_f_vect(i) = para_calc(Co_f,H_C_f,T+273);
end



%PLOT PARAMETERS V TEMP%
figure('NumberTitle', 'off', 'Name', 'GAB Model Constants v Temperature')
subplot(1,3,1)
hold on
plot(temp_profile,M0_s_vect,'-r');
plot(temp_profile,M0_f_vect,'-b');
%format
title('Monolayer, M0')
xlabel('Temperature [C]');
ylabel('M0 [ ]');
legend('Semolina','Farina');

subplot(1,3,2)
hold on
plot(temp_profile,K_s_vect,'-r');
plot(temp_profile,K_f_vect,'-b');
%format
title('K Constant')
xlabel('Temperature [C]');
ylabel('K [ ]');
legend('Semolina','Farina');

subplot(1,3,3)
hold on
plot(temp_profile,C_s_vect,'r-');
plot(temp_profile,C_f_vect,'b-');
%format
title('C Constant')
xlabel('Temperature [C]');
ylabel('C [ ]');
legend('Semolina','Farina');

%initialize the water activity mesh
num_steps = 100;
aw_mesh = zeros(1,num_steps+1);
step_size = 1/num_steps;

%populate the mesh
for i = 1:1:length(aw_mesh)-1
   aw_mesh(i+1) = aw_mesh(i) + step_size;
end

%initialize the temperatures to be used
temp_mesh = [20 35 50];

X_s = zeros(length(temp_mesh),length(aw_mesh));

%CALCULATE THE MOISTURE CONTENT AS A FUNCTION OF EATER ACTIVITY USING THE
%GAB MODEL

%Semoline, 20 C GAB MODEL
for i = 1:1:length(aw_mesh)
    
    temp_ind = 1;
    temp = temp_mesh(temp_ind) + 273;
    M0 = para_calc(M0o_s,H_M0_s,temp);
    K = para_calc(Ko_s,H_K_s,temp);
    C = para_calc(Co_s,H_C_s,temp);
    aw = aw_mesh(i);
    X_s(temp_ind,i) = (M0*K*C*aw)/((1-K*aw)*(1+(C-1)*K*aw));
    
end

%Semolina, 35 C GAB MODEL
for i = 1:1:length(aw_mesh)
    
    temp_ind = 2;
    temp = temp_mesh(temp_ind) + 273;
    M0 = para_calc(M0o_s,H_M0_s,temp);
    K = para_calc(Ko_s,H_K_s,temp);
    C = para_calc(Co_s,H_C_s,temp);
    aw = aw_mesh(i);
    X_s(temp_ind,i) = (M0*K*C*aw)/((1-K*aw)*(1+(C-1)*K*aw));
    
end

%Semoline, 50 C GAB MODEL
for i = 1:1:length(aw_mesh)
    
    temp_ind = 3;
    temp = temp_mesh(temp_ind) + 273;
    M0 = para_calc(M0o_s,H_M0_s,temp);
    K = para_calc(Ko_s,H_K_s,temp);
    C = para_calc(Co_s,H_C_s,temp);
    aw = aw_mesh(i);
    X_s(temp_ind,i) = (M0*K*C*aw)/((1-K*aw)*(1+(C-1)*K*aw));
    
end

%intialize the farina moisture content vector
X_f = zeros(length(temp_mesh),length(aw_mesh));

%FARINA, 20 C GAB MODEL
for i = 1:1:length(aw_mesh)
    
    temp_ind = 1;
    temp = temp_mesh(temp_ind) + 273;
    M0 = para_calc(M0o_f,H_M0_f,temp);
    K = para_calc(Ko_f,H_K_f,temp);
    C = para_calc(Co_f,H_C_f,temp);
    aw = aw_mesh(i);
    X_f(temp_ind,i) = (M0*K*C*aw)/((1-K*aw)*(1+(C-1)*K*aw));
    
end

%FARINA, 35 C GAB MODEL
for i = 1:1:length(aw_mesh)
    
    temp_ind = 2;
    temp = temp_mesh(temp_ind) + 273;
    M0 = para_calc(M0o_f,H_M0_f,temp);
    K = para_calc(Ko_f,H_K_f,temp);
    C = para_calc(Co_f,H_C_f,temp);
    aw = aw_mesh(i);
    X_f(temp_ind,i) = (M0*K*C*aw)/((1-K*aw)*(1+(C-1)*K*aw));
    
end

%FARINA, 50 C GAB MODEL
for i = 1:1:length(aw_mesh)
    
    temp_ind = 3;
    temp = temp_mesh(temp_ind) + 273;
    M0 = para_calc(M0o_f,H_M0_f,temp);
    K = para_calc(Ko_f,H_K_f,temp);
    C = para_calc(Co_f,H_C_f,temp);
    aw = aw_mesh(i);
    X_f(temp_ind,i) = (M0*K*C*aw)/((1-K*aw)*(1+(C-1)*K*aw));
    
end

% PLOT THE DATA THAT WAS JUST CALCUALTED. MOISTURE CONTENT AS A FUNCTION OF
% WATER ACTIVITY
figure('NumberTitle', 'off', 'Name', 'Mositure Isotherms for Product (FITTED DATA)')
subplot(2,1,1)
hold on
plot(aw_mesh,X_s(1,:));
plot(aw_mesh,X_s(2,:));
plot(aw_mesh,X_s(3,:));
%plot(aw_mesh,X_s(4,:));
legend('20 C', '35 C', '50 C','location','northwest');
title('Moisture Content v Water Activity (Semolina)');
xlabel('Water Activity, aw')
ylabel('Moisture Content');

subplot(2,1,2)
hold on
plot(aw_mesh,X_f(1,:));
plot(aw_mesh,X_f(2,:));
plot(aw_mesh,X_f(3,:));
%plot(aw_mesh,X_f(4,:));
legend('20 C', '35 C', '50 C','location','northwest');
title('Moisture Content v Water Activity (Ferina)');
xlabel('Water Activity, aw')
ylabel('Moisture Content');

%Define the temperature mesh.
temp_mesh = [20 35 50 60];

X_s_emp = zeros(length(temp_mesh),length(aw_mesh));

%CALCUALTE THE GAB MODEL EQUATION AGAIN BUT USING THE EMPIRCALLY DEFINED
%CONSTANTS
for i = 1:1:length(aw_mesh)
    
    temp_ind = 1;
    temp = temp_mesh(temp_ind) + 273;
    M0 = 11.8;
    K = 0.65;
    C = 4.21;
    aw = aw_mesh(i);
    X_s_emp(temp_ind,i) = (M0*K*C*aw)/((1-K*aw)*(1+(C-1)*K*aw));
    
end

for i = 1:1:length(aw_mesh)
    
    temp_ind = 2;
    temp = temp_mesh(temp_ind) + 273;
    M0 = 6.45;
    K = 0.71;
    C = 8.77;
    aw = aw_mesh(i);
    X_s_emp(temp_ind,i) = (M0*K*C*aw)/((1-K*aw)*(1+(C-1)*K*aw));
    
end
for i = 1:1:length(aw_mesh)
    
    temp_ind = 3;
    temp = temp_mesh(temp_ind) + 273;
    M0 = 5.92;
    K = 0.72;
    C = 10.01;
    aw = aw_mesh(i);
    X_s_emp(temp_ind,i) = (M0*K*C*aw)/((1-K*aw)*(1+(C-1)*K*aw));
    
end


X_f_emp = zeros(length(temp_mesh),length(aw_mesh));

for i = 1:1:length(aw_mesh)
    
    temp_ind = 1;
    temp = temp_mesh(temp_ind) + 273;
    M0 = 9.16;
    K = 0.65;
    C = 14.06;
    aw = aw_mesh(i);
    X_f_emp(temp_ind,i) = (M0*K*C*aw)/((1-K*aw)*(1+(C-1)*K*aw));
    
end

for i = 1:1:length(aw_mesh)
    
    temp_ind = 2;
    temp = temp_mesh(temp_ind) + 273;
    M0 = 6.01;
    K = 0.69;
    C = 31.95;
    aw = aw_mesh(i);
    X_f_emp(temp_ind,i) = (M0*K*C*aw)/((1-K*aw)*(1+(C-1)*K*aw));
    
end

for i = 1:1:length(aw_mesh)
    
    temp_ind = 3;
    temp = temp_mesh(temp_ind) + 273;
    M0 = 5.23;
    K = 0.71;
    C = 37.04;
    aw = aw_mesh(i);
    X_f_emp(temp_ind,i) = (M0*K*C*aw)/((1-K*aw)*(1+(C-1)*K*aw));
    
end

%PLOT THE GAB MODEL EQUATION DATA AS CALCULATED USING THE EXPERIMETNALLY
%DETERMINED DATA
figure('NumberTitle', 'off', 'Name', 'Mositure Isotherms for Product (ACTUAL DATA)')
subplot(2,1,1)
hold on
plot(aw_mesh,X_s_emp(1,:));
plot(aw_mesh,X_s_emp(2,:));
plot(aw_mesh,X_s_emp(3,:));
%plot(aw_mesh,X_s_emp(4,:));
legend('20 C', '35 C', '50 C','location','northwest');
title('Moisture Content v Water Activity (Semolina)');
xlabel('Water Activity, aw')
ylabel('Moisture Content');

subplot(2,1,2)
hold on
plot(aw_mesh,X_f_emp(1,:));
plot(aw_mesh,X_f_emp(2,:));
plot(aw_mesh,X_f_emp(3,:));
%plot(aw_mesh,X_f_emp(4,:));
legend('20 C', '35 C', '50 C','location','northwest');
title('Moisture Content v Water Activity (Ferina)');
xlabel('Water Activity, aw')
ylabel('Moisture Content');

%Eb calculation for Semolina
%0.1 wb = 0.1111111 db
%0.6 wb = 1.5 db

moist_vect = linspace(0.11111,1.5,100);
moist_vect_wb = zeros(1,length(moist_vect));

for i = 1:1:length(moist_vect)
    moist_vect_wb(i) = moist_vect(i)/(1 + moist_vect(i));
end

Eb_vect_S = zeros(1,length(moist_vect));
Eb_vect_F = zeros(1,length(moist_vect));

%Define 3 temperatures
T1 = 20 + 273;
T2 = 35 + 273;
T3 = 50 + 273;

%CALCULATE THE BINDING ENERGY FOR THE SEMOLINA
for i = 1:1:length(moist_vect)
    M = moist_vect(i);
    aw1 = aw_calc(M,para_calc(M0o_s,H_M0_s,T1),para_calc(Co_s,H_C_s,T1),para_calc(Ko_s,H_K_s,T1));
    aw2 = aw_calc(M,para_calc(M0o_s,H_M0_s,T2),para_calc(Co_s,H_C_s,T2),para_calc(Ko_s,H_K_s,T2));
    aw3 = aw_calc(M,para_calc(M0o_s,H_M0_s,T3),para_calc(Co_s,H_C_s,T3),para_calc(Ko_s,H_K_s,T3));
    Eb1 = (log(aw2/aw1)*R)/(1/T1 - 1/T2);
    Eb2 = (log(aw3/aw1)*R)/(1/T1 - 1/T3);
    Eb3 = (log(aw3/aw2)*R)/(1/T2 - 1/T3);
    Eb_avg = (Eb1 + Eb2 + Eb3)/3;
    Eb_vect_S(i) = Eb_avg;
end

%CALCULATE THE BINDING ENERGY FOR THE FARINA
for i = 1:1:length(moist_vect)
    M = moist_vect(i);
    aw1 = aw_calc(M,para_calc(M0o_f,H_M0_f,T1),para_calc(Co_f,H_C_f,T1),para_calc(Ko_f,H_K_f,T1));
    aw2 = aw_calc(M,para_calc(M0o_f,H_M0_f,T2),para_calc(Co_f,H_C_f,T2),para_calc(Ko_f,H_K_f,T2));
    aw3 = aw_calc(M,para_calc(M0o_f,H_M0_f,T3),para_calc(Co_f,H_C_f,T3),para_calc(Ko_f,H_K_f,T3));
    Eb1 = (log(aw2/aw1)*R)/(1/T1 - 1/T2);
    Eb2 = (log(aw3/aw1)*R)/(1/T1 - 1/T3);
    Eb3 = (log(aw3/aw2)*R)/(1/T2 - 1/T3);
    Eb_avg = (Eb1 + Eb2 + Eb3)/3;
    Eb_vect_f(i) = Eb_avg;
end


%PLOT THE BINDING ENERGY DATA AS A FUNCTION OF MOISTURE CONTENT
figure('NumberTitle', 'off', 'Name', 'Binding Energy v Moisture Content')
hold on
plot(moist_vect_wb,Eb_vect_S,'r-');
plot(moist_vect_wb,Eb_vect_f,'b-');
title('Binding Energy, J/kg');
xlabel('Moisture Content (Wet Basis)')
ylabel('Binding Energy [J/g]')
legend('Semolina','Farina');

%CALCULATE THE EFFECTIVE DIFFUSIVITY USING THE FOLLOWING EQUATION FOR SEMOLINA.:
%Deff = Do*exp(-Ea/RT)*(K*exp(-Eb/RT)/(1 + K*exp(-Eb/RT))
%Define constants
temps = [-39.62, 1.03, 41.64, 82.70];
K = 1032.6;
Ea = 5.2*4.184;
Do = 7e-8;
Deff_vect_S = zeros(length(temps),length(Eb_vect_S));

%loop through each moiture content and temperature to calcualte the
%diffusivity,.
for j = 1:1:length(temps)
    T = temps(j) + 273;;
    for i = 1:1:length(Deff_vect_S(1,:))
       Eb = Eb_vect_S(i);
       Deff_vect_S(j,i) = Do*exp(-1*Ea/(R*T))*(K*exp(-1*Eb/(R*T))/(1 + K*exp(-1*Eb/(R*T))));

    end   
end

%CALCULATE THE EFFECTIVE DIFFUSIVITY USING THE FOLLOWING EQUATION FOR FARINA.:
%Deff = Do*exp(-Ea/RT)*(K*exp(-Eb/RT)/(1 + K*exp(-Eb/RT))
%Define constants
K = 1032.6;
Ea = 5.2*4.184;
Do = 7e-8;
Deff_vect_f = zeros(length(temps),length(Eb_vect_S));

for j = 1:1:length(T_K)
    T = temps(j) + 273;
    for i = 1:1:length(Deff_vect_f(1,:))
       Eb = Eb_vect_f(i);
       Deff_vect_f(j,i) = Do*exp(-1*Ea/(R*T))*(K*exp(-1*Eb/(R*T))/(1 + K*exp(-1*Eb/(R*T))));

    end   
end

%PLOT THE DATA
figure('NumberTitle', 'off', 'Name', 'Diffusion Coefficient (Semolina)')
hold on
plot(moist_vect_wb,Deff_vect_S(1,:),'-k');
plot(moist_vect_wb,Deff_vect_S(2,:),'-r');
plot(moist_vect_wb,Deff_vect_S(3,:),'-g');
plot(moist_vect_wb,Deff_vect_S(4,:));
title('Diffusion Coefficient (Semolina)');
xlabel('Moisture Content (Wet Basis)')
ylabel('Diffusion Coefficient')
legend('-39.62 C', '1.03 C', '41.64 C', '82.70 C', 'location', 'southeast');

figure('NumberTitle', 'off', 'Name', 'Diffusion Coefficient (Farina)')
hold on
plot(moist_vect_wb,Deff_vect_f(1,:),'-k');
plot(moist_vect_wb,Deff_vect_f(2,:),'-r');
plot(moist_vect_wb,Deff_vect_f(3,:),'-g');
plot(moist_vect_wb,Deff_vect_f(4,:));
title('Diffusion Coefficient (Farina)');
xlabel('Moisture Content (Wet Basis)')
ylabel('Diffusion Coefficient')
legend('-39.62 C', '1.03 C', '41.64 C', '82.70 C', 'location', 'southeast');

Tg_soy = 410; %Kelvin
Tg_w = 134; %Kelvin

moist_vect_wet = linspace(0.05,0.4,100);
moist_vect_dry = zeros(1,length(moist_vect_wet));

for i = 1:1:length(moist_vect_wet)
   
    moist_vect_dry(i) = moist_vect_wet(i)/(1 - moist_vect_wet(i));
    
end

Tg_vect = zeros(1,length(moist_vect_wet));

for i = 1:1:length(moist_vect_wet)
    M = moist_vect_wet(i);
    Tg_vect(i) = 1/(M/Tg_w + (1-M)/Tg_soy);
end

Tg_vect = Tg_vect - 273; % Convert to celcius

T_dry_vect = ones(1,length(moist_vect_wet));
T_dry_vect = T_dry_vect.*Tg_vect(length(Tg_vect)) + 50;
drying_temps_s = [];
drying_temps_s(1,1) = Tg_vect(length(Tg_vect)) + 50;
drying_temps_s(2,1) = moist_vect_wet(length(moist_vect_wet));

for i = length(moist_vect_wet):-1:2;
    
    T_min = Tg_vect(i) + 10;
    T_dry = T_dry_vect(i);
    
    if (T_dry <= T_min)
        T_dry = Tg_vect(i) + 50;
        drying_temps_s(1,length(drying_temps_s(1,:))+1) = T_dry;
        drying_temps_s(2,length(drying_temps_s(2,:))) = moist_vect_wet(i-1);
    end
    
    T_dry_vect(i-1) = T_dry;
end

%disp(T_dry_vect)

figure('NumberTitle', 'off', 'Name', 'Glass Transition Temperature')
hold on
plot(moist_vect_wet,Tg_vect + 10,'b--');
plot(moist_vect_wet,Tg_vect + 50,'k--');
plot(moist_vect_wet,Tg_vect,'g-');
plot(moist_vect_wet,T_dry_vect,'r');
plot(drying_temps_s(2,:),drying_temps_s(1,:),'ob');
xlabel('Moisture Content (Wet Basis)');
ylabel('Tg [C]');
title('Glass Transition Temp. v Moisture Content');
legend('Minimum temperature','Maximum Temperature','Tg [C]','Drying Temperature');

%Define constants
H_M0_s = 2.195e4;
H_C_s =-6.539e4;
H_K_s =-2.899e3;
H_M0_f = 1.295e4;
H_C_f =-4.030e4;
H_K_f = -2.280e3;
M0o_s = 0.0014;
Co_s = 1.299e12;
Ko_s = 2.156;
M0o_f = 0.0424;
Co_f = 2.011e8;
Ko_f = 1.665;



delM = moist_vect_wb(2) - moist_vect_wb(1);
Deff_vect = [6.922e-7, 6.933e-7, 6.942e-7, 6.949e-7];
mass = 2000; %g
del_t = 0;
rho_w = 1000;
res_time = zeros(1,length(temps));

figure('NumberTitle', 'off', 'Name', 'Drying Time')
hold on
for i = 1:1:length(drying_temps_s)
    
    r = 5/1000; %m
    Deff = Deff_vect(i);
    Mo = drying_temps_s(2,i);
    
    if i == 4
        Mf = 0.1;
    else
        Mf = drying_temps_s(2,i+1);
    end
    
    Mi = Mo;
    
    while Mi> Mf
        
        Mi_d = Mi/(1 - Mi);
        del_t = del_t + ((4*r^2)/((pi()*Deff))*log(8*Mi_d/(Mi_d + delM)));
        %r = (r^3 - (3*delM*mass*rho_w)/(4*pi()))^(1/3);
        Mi = Mi - delM;
        plot(Mi,del_t/60/60,'ob');
        
    end
    
    res_time(i) = del_t;
end

xlabel('Moisture Content (Wet Basis)');
ylabel('Time (sec)');
title('Cumulative Drying Time');

res_time = res_time./60./60;


fprintf('A %d stage dryer using the following conditions: \n',length(drying_temps_s(1,:)));
fprintf('-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-\n');

for i = 1:1:length(drying_temps_s(1,:))
    M = drying_temps_s(2,i);
    temp = drying_temps_s(1,i) + 273;
    
    M0 = para_calc(M0o_s,H_M0_s,temp);
    C = para_calc(Co_s,H_C_s,temp);
    K = para_calc(Ko_s,H_K_s,temp);
    
    RH = aw_calc(M,M0,C,K)*100;
    
    fprintf('Stage %d: Air at %0.2f C and %0.2f%% RH for %0.2f hrs \n',i,temp-273,RH,res_time(i));
    
end

fprintf('-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-\n');

width = 8 * 0.3048; %ft to meters
belt_speed = 160; % meters per hour
belt_lengths = zeros(1,length(res_time));

for i = 1:1:length(res_time)
    belt_lengths(i) = res_time(i) * belt_speed;
end

for i = 1:1:length(belt_lengths)
    fprintf('Stage %d: Belt length is %0.2f meters\n',i,belt_lengths(i))
end

fprintf('-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-\n');


cost_conv = 0;
%calculate the cost of the conveyor belt
for i = 1:1:length(belt_lengths)
    cost_conv = cost_conv + 6*10^4*(belt_lengths(i)/60)^0.85;
end

%Calculate the gas flow rate
h = 100; %W/m^2-K
cp = 1; %J/kg-K
mu = 1.825e-5; %poise
k = 10*mean(Deff_vect); %kg/m^2
rho = 1.184; %density of something
v_air = 10; %m/s

Jd = k/v_air*(mu/(rho*mean(Deff_vect)))^(2/3);
Jh = Jd;

G = h/(cp*Jh)*(cp*mu/K)^(2/3); %kg/hr

T1 = 42;
T2 = 82.22;
T3 = 123.55;

%Caclulate the cost of heating steam
delH = 15; %kj/kg
steam_heat_1 = G*cp*(T1 - 25)/(delH);
steam_heat_2 = G*cp*(T2 - 25)/(delH);
steam_heat_3 = G*cp*(T3 - 25)/(delH);
steam_tot = steam_heat_1 + steam_heat_2 + steam_heat_3;
steam_tot * 2.2; %lbs/hr

cost_steam = (5/1000)*steam_tot;

pump_pressures = zeros(1,length(belt_lengths));
rho_f = 1.12; %kg/m^3
H = 0.1; %m
x = 0.005; %m
epsilon = 0.05; %voidage

for i = 1:1:length(pump_pressures)
    
   area = belt_lengths(i)*width; 
   pump_pressures(i) = H*1.75*rho_f*(v_air)^2*(1-epsilon)/(x*epsilon^2);
   
end

p1 = 101.3; %kpa
effic = 0.80; %efficieny
elect_price = 0.04; %$/kwhr
pump_powers = zeros(1,length(pump_pressures));
pump_work = zeros(1,length(pump_pressures));
pump_costs = zeros(1,length(pump_work));

for i = 1:1:length(pump_powers)
    
   pump_work(i) = (p1/rho_f)*log(pump_pressures(i)/p1);
   
end

for i = 1:1:length(pump_work)
    
   pump_powers(i) = pump_work(i)*G/(effic*1000);
    
end

for i = 1:1:length(pump_costs)

    pump_costs(i) = pump_powers(i)*res_time(i)*elect_price;
    
end

cost_pump = sum(pump_costs);

%calculate the air flow rate
vol_flow = G/1.125; %m^3 air per hour
cost_fan = 10^3*(vol_flow/84.9)^0.54; % $

fprintf('The cost of the conveyor belt is: $%0.2f\n',cost_conv);
fprintf('The cost of the steam is: $%0.2f per hour\n',cost_steam);
fprintf('The cost of the fan is: $%0.2f\n',cost_fan);
fprintf('The cost to pump the air is: $%0.2f\n',cost_pump);

function const = para_calc(const0,H,T);
    R = 8.324; %J/mol-K
    const = const0*exp(H/(R*T));
end

function aw = aw_calc(M,M0,C,K);
    
    A = (M0/M) - 1;
    aw = (2 + A*C - ((2 + A*C)^2 - 4*(1-C))^0.5)/(2*K*(1-C));

end


