function main
tic
%clear commmand line and variable space
clc
clear all
q = 2 * 60 * 60; %L/hr
% GET THE GUESS VOLUME FROM USER %
guess_V = input('Enter a guess volume (L): ');
error = 100; % Assume very high error at begining to enter loop

% KEEP CHECKING THE ERROR UNTIL IT IS BELOW OUR DESIRED VALUE %
iterations = 0;
while abs(error) > 0.1 % we want to get close %
    iterations = iterations + 1;
    % CALCULATE CONDITIONS BASED OFF OF GUESS VOLUME %
    S0 = 200*guess_V; %grams
    X_f = S0*0.95*0.5; %grams
    X_0 = 0.1*X_f;
    fill_time = guess_V/q;
    

    % INITIALIZE MESHES %
    X_mesh = [];
    S_mesh = [];
    t_mesh = [];

    % POPULATE MESHES WITH INITIAL VALUES %
    X_mesh(1) = X_0;
    S_mesh(1) = S0;
    t_mesh(1) = 0;
    
    % START COUNTER AND DEFINE STEP SIZE %
    cntr = 1;
    h = 0.0001;

    % SOLVE DIFF EQs % WE ARE USING EULERS METHOD TO SOLVE EQUATIONS %
    while (X_mesh(cntr) < X_f) % Keep iterating until the maximum cell count is reached
        
         % EULERS METHOD X_n+1 = X_N + h*dXdt
         X_mesh(cntr + 1) = X_mesh(cntr) + h*dXdt(X_mesh(cntr),S_mesh(cntr),t_mesh(cntr),guess_V);
         S_mesh(cntr + 1) = S_mesh(cntr) + h*dSdt(X_mesh(cntr),S_mesh(cntr),t_mesh(cntr),guess_V);
         t_mesh(cntr + 1) = t_mesh(cntr) + h;
         cntr = cntr + 1;
         
    end
    
    %Extract fermentation time
    ferment_time = t_mesh(cntr); % Last time in the mesh
    
    %Calculate the rate based on the fill and ferment times
    calc_rate = X_f/(fill_time + ferment_time); %g/hr
    calc_rate = calc_rate*0.0022; %pounds/hr
    
    % Get the amount of error. 100 pounds per hour is desired.
    error = calc_rate - 100;
    
    % Gain variable for use in convergence algorithm
    kp = 1;
    
    % Converge towards a volume that creates a smaller error %
    if error < 0 % larger volume required
        guess_V = guess_V + abs(error)*kp;
    end
    
    if error > 0 % smaller volume required
        guess_V = guess_V - abs(error)*kp;
    end

end

%Plot the data
figure(1)
plot(t_mesh,X_mesh./1000,'-k');
hold on
plot(t_mesh,S_mesh./1000,'-b');
title('Yeast Growth and Substrate Level over Time')
xlabel('Time [hrs]');
ylabel('Growth/Substrate Level [kg]')
legend('Dry Yeast Level','Substrate Level')

time = toc;

% OUTPUT RESULTS %
fprintf('-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-\n');
fprintf('Volume: %0.2f L\n',guess_V);
fprintf('Fill Time: %0.2f hrs\n',fill_time);
fprintf('Fermentation Time: %0.3f hrs \n',ferment_time);
fprintf('Calculated Rate: %0.2f #/hr\n',calc_rate);
fprintf('Error: %0.2f #/hr\n',error);
fprintf('Full Fill , Ferment, Empy Cycle Time: %0.2f hrs\n', (fill_time + ferment_time)*2);
fprintf('Elapsed Implementation Time: %0.4f sec in %d iterations\n',time,iterations);
fprintf('-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-\n');

% FERMENTER PART B %
% POWER CONSUMPTION %
D_t = (2*(guess_V/1000)/pi) ^ (1/3); % m
visc_w = 0.001; %pa.s
shear_stress = 2500; %pa
N = 2.0; %RPS
RPM = N*60; %RPM
D_i = 0.5*D_t; %m
rho_w = 1000; %kg/m^3
Re_im = (D_i^2*N*rho_w)/visc_w; %unitless
Power_im = 3; %from figure 3.4-5 geankopolis
Power = Power_im*rho_w*N^3*D_i^5;

% EXCHANGER FLOW AND SIZE %
X = 200*0.95*0.5*guess_V; %g dry weight
D_ex = 0.05; %meters (5 cm)
q_O2 = 8 * 32.02; %mg O2/g
q_cells = 0.12*q_O2*X; %kcal/hr
Cp_w = 1; %kcal/hjr
delT = 10; %K
m_dot = q_cells/(Cp_w*delT); %kg/hr
m_dot_s = m_dot/60/60; %kg/s
k_w = 0.6; %W/m-K
Re = (4*m_dot_s)/(pi*visc_w*D_ex); %Unitless
Pr = Cp_w*visc_w/k_w; %unitless
Nu = 0.027*(Re^0.8)*(Pr^(1/3)); %Unitless
h = Nu*k_w/D_ex; %W/m^2-K
T_in_f = 30; %C
T_in_w = 15; %C
T_out_f = 30; %C
T_out_w = 25; %C
delT_in = T_in_f - T_in_w; %C
delT_out = T_out_f - T_out_w; %C
delT_lm = (delT_in - delT_out)/(log(delT_in/delT_out)); %C
length = q_cells/(pi*D_ex*h*delT_lm); %m
p_O2 = 0.21; %atm
H_O2 = 4.75E4; %atm/mol
X_O2 = p_O2/H_O2; %unitless


% BLOWER RATE %
C_star = guess_V*(1/1)*1000*(1/18.018)*X_O2*32.02*1000*(1/2363); %mg/L
CL = 3; %mg/L
OUR = q_O2*X; %Uptake

T = 30; %C
RPM = N*60; %RPM
kla = OUR/(C_star-CL); %1/hr
R_gas = 0.082; %L-atm/mol-K
q_O2_mmol = 8; %mmol/hr-g
q_O2_mol = q_O2_mmol/1000; %mol/hr-g
p_O2 = 0.21;
molar_flow_rate_air = q_O2_mol*X/p_O2; %mol/hr
Q_air = molar_flow_rate_air*R_gas*(303)/1; %L/hr
Q_air_m3 = Q_air/1000; %m^3/hr
Q_air_m3_s = Q_air_m3/60/60; %m^3/s
d_turb = 0.02; %m (5 cm) 
C_turb = pi*(d_turb/2)^2; %m^2
v_air = Q_air_m3/C_turb; %m/hr
v_air_s = v_air/60/60; %m/s
P_out = 101300; %Pa
rho_air = 1.164; %kg/m3
P_in = ((v_air_s^2)*rho_air+4*P_out)/(4-(v_air_s^2)*(rho_air/P_out)); %Pa
P_in_kpa = P_in/1000; %kpa
delP = P_in - P_out;
delP_kpa = delP/1000; %kpa
gamma = 1.41; %unitless
R = 8314; %J/mol-K
T = 303; %K
M = 28.97; %g/mol
effic = 0.8;
m_dot_air = Q_air_m3_s*rho_air; %kg/s

W_s = -1*(gamma/(gamma-1))*(R*T/M)*(((P_out/P_in)^((gamma-1)/gamma))-1); %J
Pow_blower = (W_s*m_dot_air)/(effic*1000); %kW

% OUTPUT %
fprintf('');
fprintf('Power Consumption by Impeller: %0.2f Watts\n', Power);
fprintf('Heat Exchanger Length: %0.2f meters\n', length);
fprintf('Solubility of Oxygen: %0.2f mg/L\n', C_star);
fprintf('KLa: %0.2f hr^-1\n',kla);
fprintf('Flow Rate of Air: %0.2f m^3/hr\n',Q_air_m3);
fprintf('Air Velocity: %0.2f m/s\n',v_air_s);
fprintf('Pressure Differential: %0.2f kpa\n',delP_kpa);
fprintf('Shaft Work: %0.2f J/kg\n',W_s);
fprintf('Blower Power Consumption: %0.2f kW\n',Pow_blower); 
% DEFINE DERIVATIVES %

function X_slope = dXdt(X,S,t,guess_V)
% DEFINE CONSTANTS %
ks = 0.25 * guess_V; %g/L
umax = 0.5; %1/hr
Yx_s = 0.5; %x/s

u = (umax*S)/(ks + S);

X_slope = u*X;

end

function S_slope = dSdt(X,S,t,guess_V)
% DEFINE CONSTANTS %
ks = 0.25*guess_V; %g/L
umax = 0.5; %1/hr
Yx_s = 0.5; %x/s

u = (umax*S)/(ks + S);

S_slope = -1*2*u*X;

end

end

