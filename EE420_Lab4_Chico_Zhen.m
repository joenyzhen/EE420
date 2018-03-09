%% EE 420 Lab 4 Wind Resource Assessment 
clear all
close all
% Author: Joshua Cinco and Joeny Zhen 

% Equation: (v/v0) = (H/H0)^alpha;
% Where alpha is the wind shear exponent -- Also known as the friction
% coefficient

%% Part 2
% 10 Minute Interval
% Reference Velocity At Position = 38.6 m
% Wind Speed at Second Point at Position = 40 m

v_test = csvread('WindData.csv', 56, 0);
ten_t_vec = v_test(:, 4);
day = v_test(:, 3);
hours = floor(ten_t_vec/100);
minutes = ten_t_vec - (hours*100);
time = (day-1) + (hours./24) + (minutes./(24*60));
H = 49;
H0 = 38.6;
V0 = v_test(:, 5);
V_49 = v_test(:, 9);
v_ratio = V_49./V0;
h_ratio = H./H0;
alpha = log(v_ratio)./log(h_ratio);
figure(1)
plot(time(1:(end-1)), alpha(1:(end-1)));
xlabel('Time (In days)');
ylabel('Amplitude of Wind Shear Exponent');
title('Plotting Wind Shear Exponent; Ref = 38.6m and Second Point = 49m')

%% Part 3 
H_new = 59.4;
H0_new = 38.6;
V0_new = v_test(:, 5);
V_59 = v_test(:, 13);
v_ratio_new = V_59./V0_new;
h_ratio_new = H_new./H0_new;

alpha_new = log(v_ratio_new)./log(h_ratio_new);
figure(2)
plot(time(1:(end-1)), alpha_new(1:(end-1)));
xlabel('Time (In days)');
ylabel('Amplitude of Wind Shear Exponent');
title('Plotting Wind Shear Exponent; Ref = 38.6m and Second Point = 59.4m')

%% Part 4

H_4 = 59.4;
H0_4 = 38.6;
h_ratio_4 = H_4./H0_4;
v_ratio_4 = (h_ratio_4).^(alpha);
V0_4 = v_test(:, 5);
V_4 = V0_4.*v_ratio_4;
figure(3);
plot(time(1:(end-1)), V_4(1:(end-1)));
hold on
plot(time(1:(end-1)), V_49(1:(end-1)));

Error = V_4 - V_59;

figure(4);
plot(time(1:(end-1)), Error(1:(end-1)));

Error_5_4 = V_59 - V_49;
figure(5);
plot(time(1:(end-1)), Error_5_4(1:(end-1)));

%% Part 5
H_80 = 80;
H0_80 = 38.6;
h_ratio_80 = H_80./H0_80;
v_ratio_80 = (h_ratio_80).^alpha_new;
V0_80 = v_test(:, 5);
V_80 = V0_80.*v_ratio_80;
plot(time(1:(end-1)), V_80(1:(end-1)));