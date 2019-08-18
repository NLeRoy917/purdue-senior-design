
%Name ___________________

%-------------------------------------------------------------------------
%USE THIS TEMPLATE FOR EACH SCRIPT/FUNCTION THAT YOU WRITE
%For each block of code (5-10 lines) in the model provide the following:
%       HOW: How does this block of code work?
%       WHAT: What does the block of code represent?
%       WHY: Why is this block of code included in the model?
%For more details, an example of good commenting and bad commenting has
%been uploaded to Blackboard.
% -------------------------------------------------------------------------

%This block of code (Lines 16-19) is used to clear MATLAB of all variables
%and to close all currently open windows. This block of code does not
%represent anything physically in the model, it merely sets up the program
%for running this script. This block of code is included to prevent
%variable name mixups between scripts and to clean up the workspace. 

clc;%clears the command window
close all;%closes graphical and other MATLAB windows
clear all;%removes all variables and functions from the workspace

%---------------Thermal/Physical Properties of the System -----------------
%Here place any/all thermal, physical, or kinetic properties of the system.
%For each block make sure to answer the HOW/WHAT/WHY questions above. 

%This block of code (Lines 29-30) defines the initial concentration and the
%final concentration. This block of code merely creates the variables. The
%reason this block of code is added is because our loop and mathematical
%equations need a starting point and a goal to compare to. 

Co=100;%Initial Concentration of the system
Cmax=10;%The desired maximum concentration at the end of the simulation

%This block of code (Lines 37-40) defines the conditions by which the
%reaction will occur. It defines and creates the variables to be used later
%in the code. The reason this block of coded is added is to use in the
%equation during the while loop below. 
V=250;%volume in which the reaction will occur
q=15;%the reaction quotient is defined, which dictacts which way a reaction will proceed. 
k=0.0025;%the reaction constant which dictactes the reate at which the reaction occurs. 

 

%-------------------Numerical/Analytical Calculations ---------------------
%Here place any calculations, mathematical structures, looping structures,
%etc that are necessary for the modeling activity. For block of code or 
%structure make sure to answer the HOW/WHAT/WHY questions above. 

%The block of code in lines 57-62 is a while loop that continuously updates
%the concentration in the system. The reason this was used was that a while
%loop is used because there is a condition that sets the exit from the 
%loop, as opposed to a for loop that has a defined endpoint. The
%concentration needs iteratively updated, thus a loop structure was
%appropriate, as we can constantly compare the concentration to the max allowable concentration. 

while Co>Cmax %This line of code sets the ending condition to exit the loop (comparing the concentration to the max final concentration. 
    Cnew=(sqrt(q)*sqrt(4*V*Co*k+q)-q)/(2*V*k);%This equation (Eq XX-X from Geankoplis) updates the concentration in the reactor. 
    Co=Cnew;%Updates the initial concentration for the next reaction. 
    display(Co)%Displays the updating concentration
 
end

%-------------------------- Graphical/Numerical Output --------------------
%Any output parameters should go here such as display of graphs, matrices,
%values, variables, or tables. For each block of code make sure to answer the
%HOW/WHAT/WHY questions above. 

%This block of code gives the user the final answer for the system
%(concentration). This line of code displays the concentration to the
%workspace. The reason is because the final concentration is the goal of
%the program, to find the final concentration for the reaction. 
display(Co)%Displays the final value in the program of concentration


