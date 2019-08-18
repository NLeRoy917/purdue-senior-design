clc
clear all
close all

Tg_soy = 410; %Kelvin
Tg_w = 134; %Kelvin

moist_vect_wet = linspace(0.1,0.6,100);
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

figure(1)
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

fprintf('A %d stage dryer using the following conditions: \n',length(drying_temps_s(1,:)));
fprintf('-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-\n');

for i = 1:1:length(drying_temps_s(1,:))
    M = drying_temps_s(2,i);
    temp = drying_temps_s(1,i) + 273;
    
    M0 = para_calc(M0o_s,H_M0_s,temp);
    C = para_calc(Co_s,H_C_s,temp);
    K = para_calc(Ko_s,H_K_s,temp);
    
    RH = aw_calc(M,M0,C,K)*100;
    
    fprintf('Stage %d: Air at %0.2f C and %0.2f%% RH \n',i,temp-273,RH);
    
end

fprintf('-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-\n');

function const = para_calc(const0,H,T);
    R = 8.324; %J/mol-K
    const = const0*exp(H/(R*T));
end

function aw = aw_calc(M,M0,C,K);
    
    A = (M0/M) - 1;
    aw = (2 + A*C - ((2 + A*C)^2 - 4*(1-C))^0.5)/(2*K*(1-C));

end