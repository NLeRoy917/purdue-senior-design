%Name Nathan LeRoy 09/26/2018


%CLEAR THE VARIABLE SPACE:
%This block of code (Lines 10-12) is used to clear MATLAB of all variables
%and to close all currently open windows. This block of code does not
%represent anything physically in the model, it merely sets up the program
%for running this script. This block of code is included to prevent
%variable name mixups between scripts and to clean up the workspace. 
clc %clears the command line
clear all %removes all variables from the workspace
close all %closes all windows in the program

%DEFINE THE PHYSICAL COMPOSITION OF OUR PUMPKIN PIE FILLING:
%This block of code (21-26) initializes mass fractions for each food component.
%All data on physical composition was found on the USDA food database. This
%data is necesssary as it allows us to calculate the relevant thermodynamic
%properties of our food item. Without them, no analysis on the temperature
%profile of our food over time could be done

X_w = 0.48; %mass fraction of water in filling
X_pro = 0; %mass frcation of protein in filling
X_fat = 0; %mass fraction of fat in filling
X_carb = 0.511; %mass fraction of carbohydrates in filling
X_fib = 0.011; %mass fraction of fiber in filling
X_ash = 0.0004; %mass fraction of ash in filling


%DEFINE INITIAL CONSTNATS AND CONDITIONS:
%This block of code (38-41) serves to set and define the initial condtions of
%our system. This is essential for solving our system numerically. We must
%know the initial temperature profile of our can and the tempereautre of
%the steam in order to begin transient anaylsis. Hewre we will also define
%the time and size steps in both dimensions for numerical analysis. This 
%code block will also serve to initialize our matrix that will store and
%iterate through our data.

T_pie_o = 366.48; %Initial temperature of our filling, K
T_steam = 160 + 273; %steam temperature in retort, K (121 deg C)
num_slices = 10; %nuber of slices in finite difference method
num_time_steps = 1000; %number of time steps we wish to use

R = 2.125*0.0254; %can radius in meters. Convert from 2.125 inches to meters.
Z = 6/100; %can height in meters. Convert from 6 inches to meters

t_final = 10*60*60; %two hour final time, seconds
del_t = t_final/num_time_steps; %small time step for numerical analysis, seconds
del_r = R/num_slices; %small step in space in r-direction, m


%INITIALIZE OUR DATA MATRIX:
%This block of code (59-75) serves to initialize our matrix with the initial
%conditions of our system. Currently, the matrix is a zero matrix - we want
%to populate the entire thing with the initial temperature of our food
%product. Then, we want to initialize the outermost layer with our steam
%temeprature. This will serve as the starting point for the finite
%difference method.

%To optimize our code, it is necessary to predfine the size of our
%matrices. Otherwise, a large portion of our computational power is devoted
%to extending the memory of our matrices and can lead to sub-optimal
%performance when solving the finite difference method. We can do this by
%simply defining a zeros matrix that will be populated later. 
T = zeros(num_time_steps,num_slices+1);

for i = 1:1:num_time_steps %iteraete through each time step
    for j = 1:1:length(T(1,1:num_slices)) %iteraete through each slice (except for last one)
         T(i,j) = T_pie_o; %set temperature to the initial temp.
    end
end

%Initialize the outermost column to be the steam temperature
for i = 1:1:num_time_steps %iterate through each time point
    T(i,num_slices+1) = T_steam; %set outermost slice/layer to the steam temperature.
end

    for n = num_slices:-1:1 %decrementing for loop, Start at outside slice and move inward.

        % Print statements for debugging
%         fprintf('%d\n',n); %print statement for debug
%         fprintf('Density of Pie Filling: %0.2f kg/m^3\n',rho_pie);
%         fprintf('Thermal Conductivity of Pie Filling: %0.2f W/m-K\n',k_pie);
%         fprintf('Heat Capacity of Pie Filling: %0.2f J/kg-K\n',cp_pie);
        
        temp = T(i,n)-273; %convert to celsius for thermal properties
        
        therm_constants = choi_okos(X_w,X_pro,X_fat,X_carb,X_fib,X_ash,temp);

        k_pie = therm_constants(1);
        rho_pie = therm_constants(2);
        cp_pie = therm_constants(3);

        alpha = k_pie/(rho_pie*cp_pie); %thermal diffusivity of the pie filling
        M = (del_r^2)/(alpha*del_t); %define M constant for finite difference method

        if(M<4) % stability check
            error('M is not grater than 4, Unstable Solution ahead!'); %error stops script
        end
        
        if(n > 2) %Use generic finite difference equation for all points but the center
            T(i+1,n) = (1/M)*(((2*n + 1)/(2*n))*T(i,n+1) + (M-2)*T(i,n) + (2*n - 1)/(2*n)*T(i,n-1));
        else %Use the special case center-point equation.
            T(i+1,n) = (4/M)*T(i,n+1) + (M-4)/M*T(i,n);
        end
        
            
    end
end

%FINITE DIFFERENCE METHOD AND CALCULATION:
%This block of code (88-120) is the meat and potatoes of our algorithm. This nested
%loop will iterate through each time step and each slice, calculating the
%chage in temperature in both space adn time dimensions. Theoutermost loop
%iterates through each time point. At each time, we get to the next loop
%with iterates through each slice. The finite difference equation from
%geankopolis is used to get the temeprature gradient. Notice that the inner
%loop decrements. This is becasue the largest n is our outer-most slice,
%while an n of 1 corresponds to the center.

for i = 1:1:num_time_steps-1 %incrementing for loop

%CONVERT TO CELSIUS:
%Frequently, it is more intuitive and useful to view temperature  gradients
%in units of relative temperature - like cesius or fahrenheit. This simple
%for loop and block of code (129-133) iterates through each row and column 
%converting the temperature in degress Kelvin to the temperature in degrees 
%celsius by subtracting 273 from the current temperature.

for i = 1:1:length(T(:,1));
    for j = 1:1:length(T(1,:));
        T(i,j) = T(i,j) - 273;
    end
end

%CALCULATE LOG REDUCTION IN CAN:
%This code (146-169) allows us to calcualte the log reduction in the center of
%the can over time. It should be noted that we ASSUME a very long process,
%and then find the time-point where a 13.5 log reduction is achieved at the
%center. This time-point will be used as the heating process time and the
%temperature profile will be recalcualted. Without this recalculation, the
%cooling profile would be off and not accurate, as the material wuld be
%hotter than it is supposed to be. We choose to calculate the log-reduction
%here using the F0 method given our D250 value, and z-value. The units for
%this calculation are in seconds and Fahrenheit.

D250 = 0.2*60; %seconds, assume largest in range 0.1-0.2 for maximum safety
z_value = 12; %F

node = 1; % take at center
F_0_250 = 0; % intiialize at 0
F0_vect = zeros(num_time_steps); %initialize with zeros
log_red = zeros(num_time_steps,1); %initialize with zeros

for i = 1:1:num_time_steps-1 %iterate through all without exceeding matrix dimensions
    temp = T(i,node)*(9/5) + 32; %temp in Fahrenheit
    F0_vect(i+1) = F0_vect(i) + del_t*10^((temp-250)/z_value); %calciulate the F0
    log_red(i+1) = F0_vect(i)/D250; %calcualte log reduction
end

for i = 1:1:num_time_steps %iterate through each
    if log_red(i) > 13.5 %is the log-reduction greater than 13.5?
        time_heat = (i-1)*del_t/60/60; %grab time
        fprintf('A heating time of %0.2f hours produces a %0.2f log reduction\n',time_heat,log_red(i)); %print results
        break %exit loop
    end
end














% REDO CALCULATIONS WITH NEW HEATING TIME %


























%DEFINE INITIAL CONSTNATS AND CONDITIONS AGAIN:
% We will be re-conducting the finite difference calculations (215-285) with our
% new-found heating time. If this was not done, we would have an over
% heated product which would affect our cooling and nutrient degredation
% caluclations.


R = 2.125*0.0254; %can radius in meters. Convert from 2.125 inches to meters.
Z = 6/100; %can height in meters. Convert from 6 inches to meters

t_final = time_heat*60*60; %two hour final time, seconds
del_t = t_final/num_time_steps; %small time step for numerical analysis, seconds
del_r = R/num_slices; %small step in space in r-direction, m


%INITIALIZE OUR DATA MATRIX:
%This block of code (231-242) serves to initialize our matrix with the initial
%conditions of our system. Currently, the matrix is a zero matrix - we want
%to populate the entire thing with the initial temperature of our food
%product. Then, we want to initialize the outermost layer with our steam
%temeprature. This will serve as the starting point for the finite
%difference method.

T = zeros(num_time_steps,num_slices+1);

for i = 1:1:num_time_steps %iteraete through each time step
    for j = 1:1:length(T(1,1:num_slices)) %iteraete through each slice (except for last one)
         T(i,j) = T_pie_o; %set temperature to the initial temp.
    end
end

%Initialize the outermost column to be the steam temperature
for i = 1:1:num_time_steps %iterate through each time point
    T(i,num_slices+1) = T_steam; %set outermost slice/layer to the steam temperature.
end


%FINITE DIFFERENCE METHOD AND CALCULATION:
%This block of code (255-285) is the meat and potatoes of our algorithm. This nested
%loop will iterate through each time step and each slice, calculating the
%chage in temperature in both space adn time dimensions. Theoutermost loop
%iterates through each time point. At each time, we get to the next loop
%with iterates through each slice. The finite difference equation from
%geankopolis is used to get the temeprature gradient. Notice that the inner
%loop decrements. This is becasue the largest n is our outer-most slice,
%while an n of 1 corresponds to the center.

for i = 1:1:num_time_steps-1 %incrementing for loop
    for n = num_slices:-1:1 %decrementing for loop, Start at outside slice and move inward.
         %fprintf('%d\n',n); %print statement for debug
%         fprintf('Density of Pie Filling: %0.2f kg/m^3\n',rho_pie);
%         fprintf('Thermal Conductivity of Pie Filling: %0.2f W/m-K\n',k_pie);
%         fprintf('Heat Capacity of Pie Filling: %0.2f J/kg-K\n',cp_pie);
        temp = T(i,n)-273; %covnert to celsius for choi-okos
        
        therm_constants = choi_okos(X_w,X_pro,X_fat,X_carb,X_fib,X_ash,temp);

        %Extract constants%
        k_pie = therm_constants(1); 
        rho_pie = therm_constants(2);
        cp_pie = therm_constants(3);

        alpha = k_pie/(rho_pie*cp_pie); %thermal diffusivity of the pie filling
        M = (del_r^2)/(alpha*del_t); %define M constant for finite difference method

        if(M<4)
            error('M is not grater than 4, Unstable Solution ahead!'); %error ends script
        end
        
        if(n > 2) %Use generic finite difference equation for all points but the center
            T(i+1,n) = (1/M)*(((2*n + 1)/(2*n))*T(i,n+1) + (M-2)*T(i,n) + (2*n - 1)/(2*n)*T(i,n-1));
        else %Use the special case center-point equation.
            T(i+1,n) = (4/M)*T(i,n+1) + (M-4)/M*T(i,n);
        end
        
            
    end
end

%CONVERT TO CELSIUS:
%Frequently, it is more intuitive and useful to view temperature  gradients
%in units of relative temperature - like cesius or fahrenheit. This simple
%for loop and block of code (294-298) iterates through each row and column 
%converting the temperature in degress Kelvin to the temperature in degrees 
%celsius by subtracting 273 from the current temperature.

for i = 1:1:length(T(:,1));
    for j = 1:1:length(T(1,:));
        T(i,j) = T(i,j) - 273;
    end
end

%DEFINE AND INITIALIZE A MATRIX TO STORE THE COOLING TEMPERATURE:
%After we complete and model our heating process, we need to complete and
%model the cooling process. Here we imagine our cans are being submerged in
%a bath of water with a constant temperature of 55 degrees C. Our cans
%leave the retort with a temperature profile identical to the one they had
%at the end of the heating process. Thus, our cooling temeprature matrix is
%identical to our heating matrix, with the exception of the outermost layer
%temperature being 55 degrees C instead of 121 degrees C.
%Lines (311-320)


T_cool = zeros(num_time_steps,num_slices+1);
T_cool(:,num_slices+1) = 12.7;
T_cool(1,:) = T(num_time_steps,:);

%convert T_cool to kelvin for analysis
for i = 1:1:length(T_cool(:,1));
    for j = 1:1:length(T_cool(1,:));
        T_cool(i,j) = T_cool(i,j) + 273;
    end
end

%FINITE DIFFERENCE METHOD AND CALCULATION:
%We must now run the fnite difference model again with our new outer
%temperature. For simplicity, the exact same assumptions are being made and
%this, we can use the same equations. Our cans are now just being emersed
%in a cool nath of water. We run this cooling process for the same time as
%our heating process. We will then extract the time the average temperature
%in the can is 100F
%Lines(331-377).

for i = 1:1:num_time_steps-1 %incrementing for loop
    for n = num_slices:-1:1 %decrementing for loop, Start at outside slice and move inward.
         %fprintf('%d\n',n); %print statement for debug
%         fprintf('Density of Pie Filling: %0.2f kg/m^3\n',rho_pie);
%         fprintf('Thermal Conductivity of Pie Filling: %0.2f W/m-K\n',k_pie);
%         fprintf('Heat Capacity of Pie Filling: %0.2f J/kg-K\n',cp_pie);
        
        temp = T_cool(i,n)-273;
        
        therm_constants = choi_okos(X_w,X_pro,X_fat,X_carb,X_fib,X_ash,temp);

        k_pie = therm_constants(1);
        rho_pie = therm_constants(2);
        cp_pie = therm_constants(3);

        alpha = k_pie/(rho_pie*cp_pie); %thermal diffusivity of the pie filling
        M = (del_r^2)/(alpha*del_t); %define M constant for finite difference method

        if(M<4)
            error('M is not grater than 4, Unstable Solution ahead!');
        end
        
        
        if(n > 2) %Use generic finite difference equation for all points but the center
            T_cool(i+1,n) = (1/M)*(((2*n + 1)/(2*n))*T_cool(i,n+1) + (M-2)*T_cool(i,n) + (2*n - 1)/(2*n)*T_cool(i,n-1));
        else %Use the special case center-point equation.
            T_cool(i+1,n) = (4/M)*T_cool(i,n+1) + (M-4)/M*T_cool(i,n);
        end
        
%         if(n == num_slices) %Use generic finite difference equation for all points but the center
%             
%         elseif(n>2 && < num_slices %Use the special case center-point equation.
%             T(i+1,n) = (1/M)*(((2*n + 1)/(2*n))*T(i,n+1) + (M-2)*T(i,n) + (2*n - 1)/(2*n)*T(i,n-1));
%         else
%             T(i+1,n) = (4/M)*T(i,n+1) + (M-4)/M*T(i,n);
%         end
        
        
     end
end

%convert T_cool to celsius for analysis
for i = 1:1:length(T_cool(:,1));
    for j = 1:1:length(T_cool(1,:));
        T_cool(i,j) = T_cool(i,j) - 273;
    end
end

%Concatenate our heating matrix with our cooling matrix. This allows us to
%manipulate and calcualte data for the entire process instead of using two
%separate matrices.
T_process = [T;T_cool]; %concatenate on top of each other.
avg_temp = zeros(2*num_time_steps,1);

%CALCUALTE AVERAGE TEMPERATURE:
%This block of code(391-408) will calcualte the average temperature in our can for
%each time-point. This is done by weighting the temperature in each slice
%by its relative volume in the can as awhole, and then summing these
%temperatures up. It is the discrete version of the average value theorem.
%it also looks for and fidns the time the average temperature is the
%required exit temperature in our can.

for i=1:1:2*num_time_steps %run each time step
    for j = num_slices:-1:1 %decrement from edge to cetner
        avg_temp(i) = avg_temp(i) + ((1/R)^2)*T_process(i,j)*((j*del_r)^2 - ((j-1)*del_r)^2); %calcuatle avg temp
    end
end

j=0; %create counter

for i = num_time_steps:1:2*num_time_steps
    if avg_temp(i) < 37.77 %search for an average temp of 100 F
        time_cool = (j)*del_t/60/60;
        fprintf('A cooling time of %0.2f hours leaves contents of can at an average temperature of %0.2f F\n',time_cool,avg_temp(i)*9/5 + 32); %print temp in fahrenheit
        break %exit loop
    end
    j = j + 1; %increment counter
end

%CALCULATE LOG REDUCTION IN CAN:
%This code (421-433) allows us to calcualte the log reduction in the center of
%the can over time. It should be noted that we ASSUME a very long process,
%and then find the time-point where a 13.5 log reduction is achieved at the
%center. This time-point will be used as the heating process time and the
%temperature profile will be recalcualted. Without this recalculation, the
%cooling profile would be off and not accurate, as the material wuld be
%hotter than it is supposed to be. We choose to calculate the log-reduction
%here using the F0 method given our D250 value, and z-value. The units for
%this calculation are in seconds and Fahrenheit.

D250 = 0.2*60; %seconds, assume largest in range 0.1-0.2 for maximum safety
z_value = 12; %F

node = 1; % take at center
F_0_250 = 0;
F0_vect = zeros(2*num_time_steps);
log_red = zeros(2*num_time_steps,1);

for i = 1:1:2*num_time_steps-1
    temp = T_process(i,node)*(9/5) + 32; %temp in Fahrenheit
    F0_vect(i+1) = F0_vect(i) + del_t*10^((temp-250)/z_value);
    log_red(i+1) = F0_vect(i)/D250;
end

%CALCULATE VITAMIN B1 LOG REDUCTION IN CAN:
%This code (441-454) allows us to calcualte how much of vitamin B1 has been
%destroyed over time. The log-reduction specifically. This is done for each
%slice at each time point. The process is identical to the log-reduction of
%microorganisms in our can over time.

D250_v = 246.9*60; %seconds, assume largest in range 0.1-0.2 for maximum safety
z_value_v = 49; %F

node = num_slices; % take at the edge
F0v_vect = zeros(2*num_time_steps,num_slices); %initialize vector
log_red_v = zeros(2*num_time_steps,num_slices); %initialize vector

for i = 1:1:2*num_time_steps-1 %iterate through entire matrix time points
    for j = num_slices:-1:1 %iterate from outermost slice to center
        temp = T_process(i,j)*9/5 + 32; %temperature in fahrenheit
        F0v_vect(i+1,j) = F0v_vect(i,j) + del_t*10^((temp - 250)/z_value_v);
        log_red_v(i+1,j) = F0v_vect(i,j)/D250_v; %calculate log reduction
    end
end

%CALCUALTE AVERAGE LOG REDUCTION of B1:
%This block of code(391-408) will calcualte the average log-reduction in our can for
%each time-point. This is done by weighting the log-red in each slice
%by its relative volume in the can as awhole, and then summing these
%temperatures up. It is the discrete version of the average value theorem..

avg_log_red_v = zeros(2*num_time_steps,1); %intiialize vector

for i=1:1:2*num_time_steps % iterate through each time-point
    for j = num_slices:-1:1 %iterate from outermost slice to center
        avg_log_red_v(i) = avg_log_red_v(i) + ((1/R)^2)*log_red_v(i,j)*((j*del_r)^2 - ((j-1)*del_r)^2); %calcualte average
    end
end

%convert our log reduction to a percentage reduction
avg_log_red_v_perc = ones(length(avg_log_red_v),1);
for i = 1:1:length(avg_log_red_v)
    avg_log_red_v_perc(i) = (1/(10^avg_log_red_v(i)))*100;
end


%CALCULATE THE VITAMIN C REDUCTION IN THE CAN
%identical to VITAMIN B1 reduction, but use different comments.
D250_c = 1.12*24*60*60; %seconds, assume largest in range 0.1-0.2 for maximum safety
z_value_c = 52; %F

node = num_slices; % take at the edge
F0c_vect = zeros(2*num_time_steps,num_slices);
log_red_c = zeros(2*num_time_steps,num_slices);

for i = 1:1:2*num_time_steps-1
    for j = num_slices:-1:1
        temp = T_process(i,j)*9/5 + 32; %temperature in fahrenheit
        F0c_vect(i+1,j) = F0c_vect(i,j) + del_t*10^((temp - 250)/z_value_c);
        log_red_c(i+1,j) = F0c_vect(i,j)/D250_c;
    end
end

%Initialize vector to store average value
avg_log_red_c = zeros(2*num_time_steps,1);

%calcualte the average log-reduction at each
%timepoint
for i=1:1:2*num_time_steps %iterate through each time-point
    for j = num_slices:-1:1 %iterate from outermost slice to center
        avg_log_red_c(i) = avg_log_red_c(i) + ((1/R)^2)*log_red_c(i,j)*((j*del_r)^2 - ((j-1)*del_r)^2); %caluclate the average log-reduction
    end
end

%convert average log-reduction to a percent reduction
avg_log_red_c_perc = ones(length(avg_log_red_v),1);

for i = 1:1:length(avg_log_red_v)
    avg_log_red_c_perc(i) = (1/(10^avg_log_red_c(i)))*100;
end

dq = zeros(num_time_steps,1); %initialize the vector to store energy chagne


%CALCULATE ENREGY REQUIREMENTS:
%Lines (521-531) Using the discrete energy change equation Q=m*cp*delT to calcualte the
%change in energy inside our cans. This is useful for calcualting the
%economics of oru process. The choi okos equation is required to calculate
%the thermal properties in our can.9

for i = 1:1:num_time_steps
    
    temp = avg_temp(i);
    therm_constants = choi_okos(X_w,X_pro,X_fat,X_carb,X_fib,X_ash,temp);
    k_pie = therm_constants(1);
    rho_pie = therm_constants(2);
    cp_pie = therm_constants(3);
    
    dq(i) = rho_pie*pi*R^2*Z*cp_pie*(temp-(T_pie_o-273));
    
end


%PLOT THE DATA:
%This section of code (541-593) serves to only present the data visually to the user.
%This is imperitive to help the user visualize what is occuring inside our
%can over time and amake education decisions on what the optimized
%paramters and steps need to be.

time_fig = [0:del_t/60/60:2*num_time_steps*del_t/60/60 - del_t/60/60];
dist = [0:del_r*1000:R*1000];

%TEMPERATURE PROFILE AT VARIOUS TIMES
figure(1)
hold on
plot(dist,T(1,1:num_slices+1),'-k');
plot(dist,T(floor((t_final/4)/del_t + 1),1:num_slices+1),'b-');
plot(dist,T(floor((t_final/2)/del_t + 1),1:num_slices+1),'g-');
plot(dist,T(floor((3*t_final/4)/del_t + 1),1:num_slices+1),'r-');
legend('t = 0','t = 25% processing time', 't = 50% processing time', 't = 75% processing time','location','northwest');    
title(['Temperature profile in can ',num2str(t_final/60/60),' hr'])
xlabel('Distance from Center (mm)');
ylabel('Temperature (deg C)');
legend('boxoff');
%LOG REDUCTION IN C.BOT OVER TIME
figure(2)
hold on
plot(time_fig,log_red);
title('C. Botulinum Destruction | Center of Can')
ylabel('log-reduction');
xlabel('Time (hrs)');
legend('C. Botulinum','location','northwest');
legend('boxoff');
%TEMPERATURE OVER TIME AT VARIOUS LOCATIONS IN OUR CAN
figure(3)
hold on
plot(time_fig,T_process(:,1),'-r');
plot(time_fig,T_process(:,num_slices),'b-');
plot(time_fig,T_process(:,floor(num_slices/2)),'g-');
plot(time_fig,avg_temp);
title('Full Process Temperature Profile');
xlabel('Time (Hr)');
ylabel('Temperature (deg C)');
legend('0 mm','r = 50mm','r = 25 mm','Average','location','northeast');
legend('boxoff');
%VITAMIN DESTRUCTION OVER TIME
figure(4)
hold on
plot(time_fig,avg_log_red_v_perc,'b-');
plot(time_fig,avg_log_red_c_perc,'r-');
ylabel('% Nutrients in Can');
ylim([0 100]);
xlabel('Time (hrs)');
title('% Nutrient in Can');
legend('Vitamin B1','Vitamin C','location','southwest');
legend('boxoff');
%ENERGY REQUIRMENTS/ABSORPTION IN CAN
figure(5)
hold on
plot(time_fig(1:num_time_steps),dq./1000);
ylabel('Energy, kJoules');
xlabel('Time (hrs)');
title('Energy Absorption per Can');

%FUNCTION | CHOI - OKOS EQUATION:
%This function serves to calculate three thermodynamic properties of a food
%component given the mass fraction of specific items in the food, and the
%operating temperature. Specifically, it returns the value of the density
%of a food item, the thermal conductivity of the food item, and the heat
%capacity of the food item.
%INPUTS: X_w - mass fraction of water
%        X_pro - mass fraction of protein
%        X_fat - mass fraction of fat
%        X_carb - mass fraction of carbohydrates
%        X_fib - mass fraction of fiber
%        X_ash - mass fraction of ash
%OUTPUTS:
%        therms - a matrix storing the values: [k, rho, cp]
%

function therms = choi_okos(X_w, X_pro, X_fat, X_carb, X_fib, X_ash,temp_operation)

rho_w = 9.9718e2 + 3.1439e-3*temp_operation - 3.7574e-3*temp_operation^2; %desnity of water, kg/m^3
rho_pro = 1.3299e3 - 5.1840e-1*temp_operation; %desnity of protein, kg/m^3
rho_fat = 9.2559e2 - 4.1757e-1*temp_operation; %density of fat, kg/m^3
rho_carb = 1.5991e3 - 3.1046e-1*temp_operation; %density of carbohydrates, kg/m^3
rho_fib = 1.3115e3 - 3.6589-1*temp_operation; %density of fiber, kg/m^3
rho_ash = 2.4238e3 - 2.8063e-1*temp_operation; %density of ash, kg/m^3

k_w = 5.7109e-1 + 1.7625e-3*temp_operation - 6.7036e-6*temp_operation^2; %therm. cond. water, W/m-K
k_pro = 1.7881e-1 + 1.1958e-3*temp_operation - 2.7178e-6*temp_operation^2; %therm. cond. protein W/m-K
k_fat = 1.8071e-1 + 2.7604e-3*temp_operation - 1.7749e-6*temp_operation^2; %th erm. cond. fat W/m-K
k_carb = 2.0141e-1 + 1.3874e-3*temp_operation - 4.3312e-6*temp_operation^2; %therm. cond. carb W/m-K
k_fib = 1.8331e-1 + 1.2497e-3*temp_operation - 3.1683e-6*temp_operation^2; %therm. cond. fiber W/m-K
k_ash = 3.2962e-1 + 1.4011e-3*temp_operation - 2.9069e-6*temp_operation^2; %therm. cond. ash W/m-K

cp_w = 4.1762e3 - 9.0864e-2*temp_operation + 5.4731e-3*temp_operation^2; %heat capacity of water, J/kg-K
cp_pro = 2.0082e3 + 1.2089*temp_operation - 1.3129e-3*temp_operation^2; %heat capacity of protein, J/kg-K
cp_fat = 1.9842e3 + 1.4733*temp_operation - 4.8008e-3*temp_operation^2; %heat capacity of fat, J/kg-K
cp_carb = 1.5488e3 + 1.9625*temp_operation - 5.9399e-3*temp_operation^2; %heat capacity of carbs, J/kg-K
cp_fib = 1.8459e3 + 1.8306*temp_operation - 4.6509e-3*temp_operation^2; %heat capacity of fiber., J/kg-K
cp_ash = 1.0926e3 + 1.8896*temp_operation - 3.6817e-3*temp_operation^2; %heat capcity of ash, J/kg-K

%calculate the thermal properites of our food according to the choi - okos
%equation
k = X_w*k_w + X_pro*k_pro + X_fat*k_fat + X_carb*k_carb + X_fib*k_fib + X_ash*k_ash; %W/m-K
rho = X_w*rho_w + X_pro*rho_pro + X_fat*rho_fat + X_carb*rho_carb + X_fib*rho_fib + X_ash*rho_ash; %kg/m^3
cp = X_w*cp_w + X_pro*cp_pro + X_fat*cp_fat + X_carb*cp_carb + X_fib*cp_fib + X_ash*cp_ash; %J/kg-K

therms = [k,rho,cp];

end
