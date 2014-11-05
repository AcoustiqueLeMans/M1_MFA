%   Program impedance_tube_fit_init.m
%
%   Copyright (C) 2014 LAUM UMR CNRS 6613 (France)
% 	   Olivier DAZEL <olivier.dazel@univ-lemans.fr>
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

%
% This program is associated to the course 
% "Mathematics for Acoustics and Acousticians"
%  given for students of the MSc Master Program in Acoustics and 
%  IMDEA students
% 

clear all
close all
clc

f=440;
c=342;
k=2*pi*f/c;

x=load('position.txt');
P_pure=load('P_pure.txt');
P_pure=P_pure(:,1)+i*P_pure(:,2);
P_background=load('P_background.txt');
P_background=P_background(:,1)+i*P_background(:,2);
P_noisy=load('P_noisy.txt');
P_noisy=P_noisy(:,1)+i*P_noisy(:,2);


figure(1)
hold on
plot(x,abs(P_pure),'Linewidth',3)
plot(x,abs(P_background),'k','Linewidth',3)
plot(x,abs(P_noisy),'m','Linewidth',3)
set(gca,'Fontsize',15)
title('Modulus of the three signals')
xlabel('Position (m)')
ylabel('Amplitude (au)')

figure(2)
hold on
plot(x,angle(P_pure),'Linewidth',3)
plot(x,angle(P_background),'m','Linewidth',3)
plot(x,angle(P_noisy),'k','Linewidth',3)
set(gca,'Fontsize',15)
title('Phase of the three signals')
xlabel('Position (m)')
ylabel('Phase (rad)')

