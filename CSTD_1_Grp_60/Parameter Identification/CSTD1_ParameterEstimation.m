% #########################################################################
% TUHH :: Institute for Control Systems :: Control Lab
% #########################################################################
% Experiment CSTD1: Identification and Control of a Torsional Plant
%
% Copyright Herbert Werner and Hamburg University of Technology, 2014
% #########################################################################
% This file is to be completed by the student.
% The completed version is to be published using 
%   publish('CSTD1_ParameterEstimation','pdf') 
% and submitted as a pdf-file at least one week prior to the scheduled date
% for the experiment
%
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% !!!  The gaps in the code are denoted by XXX !!!
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%
% HINT 1:
% if you want to find out more about a certain command, just type 
% 'help command' into the matlab window
% HINT 2:
% use section evaluations (Ctrl+Enter) to run the code within a single 
% section

%----------------------------
% v.1.2 - 30-10-2017
% by Antonio Mendez G
%----------------------------
%%
clc
clear all
close all

%% 1 Input and Output Data
%
% Generate 2 sets of input and output data using the simulink file 
% CSTD1_OpenLoopModel.slx
% Save the data as .mat files containing:
% 'Plant_input_u', 'Theta1', 'Theta2', 'Theta3'.
% Dimensions must be [Nm] and [deg], respectively.
% A sampling time Ts = 1ms is assumed.
% set1: 'estimation_file', data to be used for parameter identification
% set2: 'validation_file', data to be used for system validation

estimation_file = 'Data_E.mat';
validation_file = 'Data_V.mat';

load(estimation_file);

%% 2 Measurement Matrix
%
% Complete the gaps inside the function CSTD1_BuildMesurementMatrix.m
% This function computes the torque signal T and the measurement matrix M.
% The following varibles must be in the Workspace
% Plant_input_u : excitation signals -- Torque [Nm]
% Theta1, Theta2, Theta3: measured angular positions of 3 disks [deg]

CSTD1_BuildMeasurementMatrix

%% 3 Parameter Identification
%
% Compute the parameter vector p using M and T.
% p: the parameter vector, where p = [J1 J2 J3 d1 d2 d3 k1 k2]'

%p = (inv(M'*M)*M'*Torque);
p = pinv(M)*Torque;
J1 = p(1);
J2 = p(2);
J3 = p(3);
d1 = p(4);
d2 = p(5);
d3 = p(6);
k1 = p(7);
k2 = p(8);

disp(p);

%% 4 State Space Construction
%
% The function below constructs the state space of the torsional plant
sys = CSTD1_StateSpaceConstruction(J1,J2,J3,d1,d2,d3,k1,k2);

% This file will be required in the controller design part
save('IdentifiedSystem.mat','sys')

%% 5 Validation
% 
% The validation plot is shown, this plot contains the actual torque
% and the estimated torque of a different set of signals,
% i.e. T_validation vs M_validation*p,
% where p is the identified parameter vector

sample_span = length(Theta1)-1;

time = 0:0.001:length(Torque(1:sample_span-1))/1000;

% Load the validation set of data
load(validation_file);

% Compute the M matrix for validation purposes
CSTD1_BuildMeasurementMatrix

% Name the validation variables
T_validation = Torque(1:sample_span);
M_validation = M(1:sample_span,:);

estimated_Torque = M_validation*p;
error = (estimated_Torque - T_validation).^2; %update

subplot(3,1,1)
plot(time,estimated_Torque,'r')
hold on
plot(time,T_validation);
grid on
title('Torques')
legend('Estimated Torque (M*p)', 'Actual Torque (T)')
xlabel('Time [s]')
ylabel('Torque [Nm]')

subplot(3,1,2)
plot(time, error)
grid on
title('Squared error: (T-Mp)^2')
xlabel('Time [s]')
ylabel('Error')

subplot(3,1,3)
plot(time,estimated_Torque,'r')
hold on
plot(time,T_validation);
grid on
title('Torques (Zoom)')
legend('Estimated Torque (M*p)', 'Actual Torque (T)')
xlabel('Time [s]')
ylabel('Torque [Nm]')
set(gca,'xlim',[0 1])

%% 6 Torsional Springs Safety

% 
% The plot of the twist angle of the torsional springs
% is shown. The absolute twist value should be always 
% smaller than  20 degrees.
figure()

load(estimation_file);
sample_span = length(Theta1)-1;
time = 0:0.001:length(Theta1(1:sample_span))/1000;
suptitle('All signals must be smaller than 20 degrees')

subplot(2,2,1)
plot(time, abs(Theta1-Theta2));
grid on
ylabel('Twist [deg]')
legend('Estimation |\theta_1 - \theta_2|')

subplot(2,2,2)
plot(time, abs(Theta2-Theta3));
grid on
legend('Estimation |\theta_2 - \theta_3|')

load(validation_file);
sample_span = length(Theta1)-1;
time = 0:0.001:length(Theta1(1:sample_span))/1000;

subplot(2,2,3)
plot(time, abs(Theta1-Theta2),'r');
grid on
xlabel('Time [s]')
ylabel('Twist [deg]')
legend('Validation |\theta_1 - \theta_2|')

subplot(2,2,4)
plot(time, abs(Theta2-Theta3),'r');
grid on
xlabel('Time [s]')
legend('Validation |\theta_2 - \theta_3|')

