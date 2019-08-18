% OPTIMAL TEMPERATURE FINDING%

R_gas = 8.314; %J/mol-K

%Define C. Bot parameters%
z_Cbot = 5*(47-32)/9 + 273; %Fahrenheit to Kelvin
Ea_Cbot = 64 * 4182; %J/mol

%Define Vitamin B1 Paramters
z_b1 = 5*(49-32)/9 + 273; %Fahrenheit to Kelvin
Ea_b1 = 27 * 4282; %J/mol

%Define Vitamin C Paramters
z_C = 5*(52-32)/9 + 273; %Fahrenheit to Kelvin
Ea_C = 24 * 4182; %J/mol

lnk_Cbot = zeros(1,101);
lnk_b1 = zeros(1,101);
lnk_C = zeros(1,101);

Temp = [0 + 273:1:100+273];
disp(Temp)
for i = 1:1:101
    lnk_Cbot(i) = -1*Ea_Cbot/R_gas/Temp(i);
end

for i = 1:1:101
    lnk_b1(i) = -1*Ea_b1/R_gas/Temp(i);
end

for i = 1:1:101
    lnk_C(i) = -1*Ea_C/R_gas/Temp(i);
end

for i = 1:1:101
    Temp(i) = 1/Temp(i);
end
disp(Temp)
figure(1)
hold on
plot(Temp,-1.*lnk_Cbot,'r-');
plot(Temp,-1.*lnk_b1,'b-');
plot(Temp,-1.*lnk_C,'--');
title('ln(k) v 1/T');
xlabel('1/T')
ylabel('lnk');
legend('C. Botulinum','Vitamin B1','Vitamin C','location','east');
legend('boxoff');