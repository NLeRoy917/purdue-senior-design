    clear all;
    close all;
    clc;


    % Define Constants
    flow_0 = 1230 * 1.0515e-6; % m^3/s
    mass_frac = 0.0123;
    mass_frac_out = 0.1023;
    D = 15*1e-3; % diameter of tubes, meters
    L = 2.46; % meter
    delP_0 = 50.23; % atm
    Aw = 4.01e-4; %kg H2O/s-m^2-atm
    As = 1.23e-7; % m/s
   
    % loop constants
    filters = 7;
    sections = 1;
    tubes_per = 500;
    num_steps = 1000;
    step_size = L/1000;
    total = 0;
    
    while mass_frac < mass_frac_out
     
        %define loops variables
        mass_frac = 0.0123;
        mem_loc = 0;
        v = 1.7; % m/s
        q_in = flow_0;
        flow_rate_store = [];
        flow_rate_store(1) = flow_0;
        mass_frac_store = [];
        mass_frac_store(1) = mass_frac;
        length_store = [];
        length_store(1) = mem_loc;
        delP_store = [];
        delP_store(1) = delP_0; %atm
        i = 2;
        area = 3.14*D*step_size*filters*tubes_per; % m^2
        
        while mem_loc < L
            %calculate values
            rho = calc_density(mass_frac); % kg/m^3
            c_1 = mass_frac * rho; % kg glucose/ m^3 solution
            visc = calc_visc(mass_frac,rho); % Pa-s
            delP = calc_delP(visc,rho,step_size,v,D)/10; % atm
            delPi = calc_delPi(mass_frac,rho); % atm
            Nw = Aw*(delP - delPi); % kg H2O/m^2-s
            q_m = Nw * area * (1/1000); % m^3/s
            
            %calc final variables
            q_out = q_in - q_m; % m^3/s
            mass_frac_new = mass_frac*q_in/(mass_frac*q_in + (1-mass_frac)*q_in - (Nw*area/1000)); % kg glucose/ kg solution
            mem_loc_new = mem_loc + step_size;
            v = q_out/(3.14159*(D/2)^2);
        
            %store values
            flow_rate_store(i) = q_out;
            mass_frac_store(i) = mass_frac_new;
            length_store(i) = mem_loc_new;
            delP_store(i) = delP;
        
            % update parameters
            q_in = q_out; % update the q_in for the next loop
            mass_frac = mass_frac_new;
            mem_loc = mem_loc_new;
            i = i + 1;
        end
        
        total = total + filters;
        sections = sections + 1;
        filters = filters - 1;
        
    end
    


    
%     fprintf('Final mass fraction: %0.2f %% w/w using %d tubes\n', mass_frac*100, num_tubes);
%     fprintf('Number of filter apparatas: %.0f \n', num_filters);
    
    figure(1)
    plot(length_store./(max(length_store)),mass_frac_store.*100, 'r-')
    xlabel('Process Completion');
    ylabel('Glucose Conc. (Weight %)');
    legend('Glucose Conc.');
    
    figure(2)
    plot(length_store./(max(length_store)),flow_rate_store.*3600,'b-')
    xlabel('Process Completion');
    ylabel('Flow Rate (m^3/hr)');
    legend('Flow Rate');
    
    figure(3)
    plot(length_store./(max(length_store)), delP_store,'g-');
    xlabel('Process Completion');
    ylabel('Pressure Drop Across Membrane (atm)');
    legend('Pressure (atm)');
     
    
function rho = calc_density(mass_frac)
    % Define constants
    temp = 25; % degrees C
    rho_w = 1000; % kg/m^3
    rho_g = 1586.2; % kg/m^3
    
    %calculate desnity
    rho = rho_g*mass_frac + rho_w*(1-mass_frac);
    
end

function visc = calc_visc(mass_frac,rho)
    %define constants
    MW_glucose = 0.180156; % kg glucose / mol
    
    %calculate molarity of solution
    molarity = mass_frac/rho/MW_glucose;
    
    %calculate viscosity of glucose solution
    visc = 0.95*exp(molarity)-0.006; % Pa-s

end

function delP = calc_delP(visc,rho,step_size,v,D);
    
    delP = 2*0.079*visc^0.25*rho^0.75*step_size*v^1.75*(1/D^1.25);
    
end

function delPi = calc_delPi(mass_frac,rho)
    % Define constants
    MW_glucose = 0.180156; % kg glucose/mol
    R = 8.205e-5; % m3-atm/mol-K
    temp = 25 + 273; % Kelvin
    
    %claculate the concentration of glucose
    conc_g = mass_frac/rho/MW_glucose; % mol/m^3
    
    %calculate the osmotic pressure of the solution
    delPi = (conc_g/MW_glucose)*R*temp; % atm

end

