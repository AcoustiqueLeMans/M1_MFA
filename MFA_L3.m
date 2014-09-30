%   Program MFA_L3.m
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
% This program plots the figures of Lecture 3 of the Course 
% "Mathematics for Acoustics and Acousticians"
%  given for students of the MSc Master Program in Acoustics and 
%  IMDEA students
% 


clear all 
close all
clc

% Number of values in the sequence 
n=60;
% Vector of time values
t=linspace(1,2,n)';

% Parameter of the model to recover
alpha=1
beta=2


% Parameters relative to the noise on data
noise_max=max(abs(alpha*t+beta))
NSR=0.08;
noise=NSR*noise_max;
sigma=noise*rand(n,1);
% Generation of dat with noise
y=alpha*t+beta+noise*(rand(n,1)-0.5);

% Application of the basic LMS procedure
M=[t,ones(n,1)];
X_LMS=(M'*M)\(M'*y)
y_LMS=X_LMS(1)*t+X_LMS(2);
figure
plot(t,y,'.','Markersize',15)
set(gca,'Fontsize',25)
print -djpeg L3_LMS_data.jpg
hold on
plot(t,y_LMS,'m', 'Linewidth',3,'Color','m')
set(gca,'Fontsize',25)
xlim([0.9 2.1])
ylim([2.8 4.1])
print -djpeg L3_LMS.jpg


X1=t-mean(t);
X2=y-mean(y);
disp('Coorelation Coefficient')
R=(X1'*X2/(norm(X1)*norm(X2)))


figure









% Number of values in the sequence 
n=10;
% Vector of time values
t=linspace(1,2,n)';

% Parameters relative to the noise on data
noise_max=max(abs(alpha*t+beta))
NSR=0.02;
noise=NSR*noise_max;
sigma=noise*rand(n,1);
% Generation of dat with noise
y=alpha*t+beta+noise*(rand(n,1)-0.5);

% Modification of the 5th entry to make it being abberation
y(5)=alpha*t(5)+beta-0.3;
sigma(5)=2*0.3;

% Application of the basic LMS procedure
M=[t,ones(n,1)];
X_LMS=(M'*M)\(M'*y);
y_LMS=X_LMS(1)*t+X_LMS(2);

figure
errorbar(t,y,sigma)
set(gca,'Fontsize',25)
print -djpeg L3_LMS_data.jpg
hold on
xlim([0.95*min(t) 1.05*max(t)])
ylim([0.95*min(y) 1.05*max(y)])
print -djpeg100 L3_data.jpg
plot(t,y_LMS,'m', 'Linewidth',3)
set(gca,'Fontsize',25)
xlim([0.9 2.1])
ylim([2.8 4.1])
print -djpeg L3_LMS_basic.jpg

% Application of the basic LMS procedure
% with the cancellation of the 5th entry
t_suppress=t;
y_suppress=y;
t_suppress(5)=[];
y_suppress(5)=[];
M_suppress=[t_suppress,ones(n-1,1)];
X_LMS_suppress=(M_suppress'*M_suppress)\(M_suppress'*y_suppress);
y_LMS_suppress=X_LMS_suppress(1)*t+X_LMS_suppress(2);

% Application of the LMS with the
% control by standard deviation

w=diag(1./(sigma.^2));
M_p=sqrt(w)*M;
y_p=sqrt(w)*y;
X_LMS_sigma=(M_p'*M_p)\(M_p'*y_p)
y_LMS_sigma=X_LMS_sigma(1)*t+X_LMS_sigma(2);


plot(t,y_LMS_suppress,'r', 'Linewidth',3)
plot(t,y_LMS_sigma,'g', 'Linewidth',3)
xlim([0.9 2.1])
ylim([2.8 4.1])
print -djpeg L3_LMS_all.jpg




