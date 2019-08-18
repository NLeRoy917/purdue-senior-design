time_vect = linspace(-1,17,19);
DO = [3.3, 3.3, 2.4, 1.3, 0.3, 0.1, 0.0, 0.0, 0.3, 1.0,1.6,2.0,2.4,2.7,2.9,3.0,3.1,3.2,3.2];
C_star = 7.3 % mg/L

figure(1)
plot(time_vect,DO,'-ob')
title('Oxygen Content v Time')
xlabel('Time [sec]')
ylabel('Dissolved O2 [mg/L]')
legend()

% The Oxygen Uptake Rate (OUR) is the slope of points:
    % 2, 3 , 4 , 5

