%% ASEN 3128 - Lab 1 - Group 4 - Problem 1 - Main
%
% This is the script for problem 1 of ASEN 3128 Lab 1
% The purpose of this script is to use ode45 to simulate a set of 3
% equations and then plot the positions found from the simulation as a 3x1
% subplot
%
% Authors: Ryan Charles, Cole MacPherson, and Zachary Vanlangendonck
% Class: ASEN 3128 - Aircraft Dynamics
% File: main.m
% Date Modified: 22nd January 2021
%
%
%% Housekeeping
clc;
clear;
close all;

%% Initial Conditions

tspan = [0 3];
y_0 = 0.1 + zeros(1,3);

%% ode45

[t_out,y_out] = ode45(@(t_out,y_out) odefun(t_out,y_out),tspan,y_0);

%% Plots

figure
subplot(3,1,1)
plot(t_out,y_out(:,1),'linewidth',2);
title('x position vs time');
ylabel('x');
legend(strcat('Initial Condition = ',num2str(y_0(1))),'location','northwest');
subplot(3,1,2)
plot(t_out,y_out(:,2),'linewidth',2);
title('y position vs time');
ylabel('y');
legend(strcat('Initial Condition = ',num2str(y_0(2))),'location','northwest');
subplot(3,1,3)
plot(t_out,y_out(:,3),'linewidth',2);
title('z position vs time');
ylabel('z');
xlabel('time');
legend(strcat('Initial Condition = ',num2str(y_0(3))),'location','northwest');
