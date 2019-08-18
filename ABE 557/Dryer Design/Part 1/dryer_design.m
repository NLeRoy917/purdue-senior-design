%Clear command window, clear variable workspace, and close all figures
clc
clear all
close all

%Semolina Constants%
M0s = [11.8, 6.45, 5.92, 3.53];
Ks = [0.65, 0.71, 0.72, 0.76];
Cs = [4.21, 8.77, 10.1, 200.1];


%Farina Constants
M0f = [9.16, 6.01, 5.23, 4.76];
Kf = [0.65, 0.69, 0.71, 0.73];
%Create temperature Vector
Cf = [14.06, 31.95, 37.04, 137.0];

T_C = [20, 35, 50, 60];
T_K = T_C + 273;

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

M0_s_vect = zeros(1,length(temp_profile));
K_s_vect = zeros(1,length(temp_profile));
C_s_vect = zeros(1,length(temp_profile));
M0_f_vect = zeros(1,length(temp_profile));
K_f_vect = zeros(1,length(temp_profile));
C_f_vect = zeros(1,length(temp_profile));

for i = 1:1:length(M0_s_vect);
    T = temp_profile(i);
    M0_s_vect(i) = para_calc(M0o_s,H_M0_s,T+273);
end

for i = 1:1:length(K_s_vect);
    T = temp_profile(i);
    K_s_vect(i) = para_calc(Ko_s,H_K_s,T+273);
end

for i = 1:1:length(C_s_vect);
    T = temp_profile(i);
    C_s_vect(i) = para_calc(Co_s,H_C_s,T+273);
end

for i = 1:1:length(M0_f_vect);
    T = temp_profile(i);
    M0_f_vect(i) = para_calc(M0o_f,H_M0_f,T+273);
end

for i = 1:1:length(K_f_vect);
    T = temp_profile(i);
    K_f_vect(i) = para_calc(Ko_f,H_K_f,T+273);
end

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
title('Monolayer, M0')
xlabel('Temperature [C]');
ylabel('M0 [ ]');
legend('Semolina','Farina');

subplot(1,3,2)
hold on
plot(temp_profile,K_s_vect,'-r');
plot(temp_profile,K_f_vect,'-b');
title('K Constant')
xlabel('Temperature [C]');
ylabel('K [ ]');
legend('Semolina','Farina');

subplot(1,3,3)
hold on
plot(temp_profile,C_s_vect,'r-');
plot(temp_profile,C_f_vect,'b-');
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

for i = 1:1:length(aw_mesh)
    
    temp_ind = 1;
    temp = temp_mesh(temp_ind) + 273;
    M0 = para_calc(M0o_s,H_M0_s,temp);
    K = para_calc(Ko_s,H_K_s,temp);
    C = para_calc(Co_s,H_C_s,temp);
    aw = aw_mesh(i);
    X_s(temp_ind,i) = (M0*K*C*aw)/((1-K*aw)*(1+(C-1)*K*aw));
    
end

for i = 1:1:length(aw_mesh)
    
    temp_ind = 2;
    temp = temp_mesh(temp_ind) + 273;
    M0 = para_calc(M0o_s,H_M0_s,temp);
    K = para_calc(Ko_s,H_K_s,temp);
    C = para_calc(Co_s,H_C_s,temp);
    aw = aw_mesh(i);
    X_s(temp_ind,i) = (M0*K*C*aw)/((1-K*aw)*(1+(C-1)*K*aw));
    
end
for i = 1:1:length(aw_mesh)
    
    temp_ind = 3;
    temp = temp_mesh(temp_ind) + 273;
    M0 = para_calc(M0o_s,H_M0_s,temp);
    K = para_calc(Ko_s,H_K_s,temp);
    C = para_calc(Co_s,H_C_s,temp);
    aw = aw_mesh(i);
    X_s(temp_ind,i) = (M0*K*C*aw)/((1-K*aw)*(1+(C-1)*K*aw));
    
end
% for i = 1:1:length(aw_mesh)
%     
%     temp_ind = 4;
%     temp = temp_mesh(temp_ind) + 273;
%     M0 = para_calc(M0o_s,H_M0_s,temp);
%     K = para_calc(Ko_s,H_K_s,temp);
%     C = para_calc(Co_s,H_C_s,temp);
%     aw = aw_mesh(i);
%     X_s(temp_ind,i) = (M0*K*C*aw)/((1-K*aw)*(1+(C-1)*K*aw));
%     
% end



X_f = zeros(length(temp_mesh),length(aw_mesh));

for i = 1:1:length(aw_mesh)
    
    temp_ind = 1;
    temp = temp_mesh(temp_ind) + 273;
    M0 = para_calc(M0o_f,H_M0_f,temp);
    K = para_calc(Ko_f,H_K_f,temp);
    C = para_calc(Co_f,H_C_f,temp);
    aw = aw_mesh(i);
    X_f(temp_ind,i) = (M0*K*C*aw)/((1-K*aw)*(1+(C-1)*K*aw));
    
end

for i = 1:1:length(aw_mesh)
    
    temp_ind = 2;
    temp = temp_mesh(temp_ind) + 273;
    M0 = para_calc(M0o_f,H_M0_f,temp);
    K = para_calc(Ko_f,H_K_f,temp);
    C = para_calc(Co_f,H_C_f,temp);
    aw = aw_mesh(i);
    X_f(temp_ind,i) = (M0*K*C*aw)/((1-K*aw)*(1+(C-1)*K*aw));
    
end

for i = 1:1:length(aw_mesh)
    
    temp_ind = 3;
    temp = temp_mesh(temp_ind) + 273;
    M0 = para_calc(M0o_f,H_M0_f,temp);
    K = para_calc(Ko_f,H_K_f,temp);
    C = para_calc(Co_f,H_C_f,temp);
    aw = aw_mesh(i);
    X_f(temp_ind,i) = (M0*K*C*aw)/((1-K*aw)*(1+(C-1)*K*aw));
    
end

% for i = 1:1:length(aw_mesh)
%     
%     temp_ind = 4;
%     temp = temp_mesh(temp_ind) + 273;
%     M0 = para_calc(M0o_f,H_M0_f,temp);
%     K = para_calc(Ko_f,H_K_f,temp);
%     C = para_calc(Co_f,H_C_f,temp);
%     aw = aw_mesh(i);
%     X_f(temp_ind,i) = (M0*K*C*aw)/((1-K*aw)*(1+(C-1)*K*aw));
%     
% end

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

temp_mesh = [20 35 50 60];

X_s_emp = zeros(length(temp_mesh),length(aw_mesh));

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
% for i = 1:1:length(aw_mesh)
%     
%     temp_ind = 4;
%     temp = temp_mesh(temp_ind) + 273;
%     M0 = 3.53;
%     K = 0.76;
%     C = 200.1;
%     aw = aw_mesh(i);
%     X_s_emp(temp_ind,i) = (M0*K*C*aw)/((1-K*aw)*(1+(C-1)*K*aw));
%     
% end



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

% for i = 1:1:length(aw_mesh)
%     
%     temp_ind = 4;
%     temp = temp_mesh(temp_ind) + 273;
%     M0 = 4.76;
%     K = 0.73;
%     C = 137.0;
%     aw = aw_mesh(i);
%     X_f_emp(temp_ind,i) = (M0*K*C*aw)/((1-K*aw)*(1+(C-1)*K*aw));
%     
% end

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
moist_vect= linspace(5.5,30,50);
Eb_vect = zeros(1,length(moist_vect));

T1 = 20 + 273;
T2 = 35 + 273;
T3 = 50 + 273;

for i = 1:1:length(moist_vect)
    M = moist_vect(i);
    aw1 = aw_calc(M,para_calc(M0o_s,H_M0_s,T1),para_calc(Co_s,H_C_s,T1),para_calc(Ko_s,H_K_s,T1));
    aw2 = aw_calc(M,para_calc(M0o_s,H_M0_s,T2),para_calc(Co_s,H_C_s,T2),para_calc(Ko_s,H_K_s,T2));
    aw3 = aw_calc(M,para_calc(M0o_s,H_M0_s,T3),para_calc(Co_s,H_C_s,T3),para_calc(Ko_s,H_K_s,T3));
    Eb1 = (log(aw2/aw1)*R)/(1/T1 - 1/T2);
    Eb2 = (log(aw3/aw1)*R)/(1/T1 - 1/T3);
    Eb3 = (log(aw3/aw2)*R)/(1/T2 - 1/T3);
    Eb_avg = (Eb1 + Eb2 + Eb3)/3;
    Eb_vect(i) = Eb_avg;
end

figure('NumberTitle', 'off', 'Name', 'Binding Energy v Moisture Content (Semlina)')
plot(moist_vect,Eb_vect,'b-');
title('Binding Energy, J/kg');
xlabel('Moisture Content [g/g]')
ylabel('Binding Energy [J/g]')




function const = para_calc(const0,H,T);
    R = 8.324; %J/mol-K
    const = const0*exp(H/(R*T));
end

function aw = aw_calc(M,M0,C,K);
    
    A = (M0/M) - 1;
    aw = (2 + A*C - ((2 + A*C)^2 - 4*(1-C))^0.5)/(2*K*(1-C));

end


