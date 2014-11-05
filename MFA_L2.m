%   Program MFA_L2.m
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
% This program plots the figures of Lecture 2 of the Course 
% "Mathematics for Acoustics and Acousticians"
%  given for students of the MSc Master Program in Acoustics and 
%  IMDEA students
% 


clear all 
close all
clc


figure
t=linspace(0,2*pi,300);
plot(cos(t),sin(t),'Linewidth',24)
hold on
n=4;
delta_m=0;
line([0 cos(delta_m)],[0 sin(delta_m)],'Linewidth',24,'Color','r')
delta_m=pi;
line([0 cos(delta_m)],[0 sin(delta_m)],'Linewidth',24,'Color','r')
for m=1:n
  delta_m=-m*pi/(n+1);
  line([0 cos(delta_m)],[0 sin(delta_m)],'Linewidth',24,'Color','g')
end
for m=1:n
  delta_m=m*pi/(n+1);
  line([0 cos(delta_m)],[0 sin(delta_m)],'Linewidth',24)
end
axis equal
axis off
print -djpeg circle_even.jpg

figure
t=linspace(0,2*pi,300);
plot(cos(t),sin(t),'Linewidth',24)
hold on
n=3;
delta_m=0;
line([0 cos(delta_m)],[0 sin(delta_m)],'Linewidth',24,'Color','r')
delta_m=pi;
line([0 cos(delta_m)],[0 sin(delta_m)],'Linewidth',24,'Color','r')
for m=1:n
  delta_m=-m*pi/(n+1);
  line([0 cos(delta_m)],[0 sin(delta_m)],'Linewidth',24,'Color','g')
end
for m=1:n
  delta_m=m*pi/(n+1);
  line([0 cos(delta_m)],[0 sin(delta_m)],'Linewidth',24)
end
axis equal
axis off
print -djpeg circle_odd.jpg


% Number of defrees of freedom 
n=50;
% Number of the selected modes
modes_list=[1 4 30 n];

% Stiffness and Mass matrices
K=zeros(n,n);
M=zeros(n,n);
% Values for the stiffness of the springs
k=1;
% Values for the individual masses
m=1;
omega_0=sqrt(k/m);


% Values in the global matrices
K(1,1)=2*k;
M(1,1)=m;

for ii=2:n
    K(ii,ii-1) = -k;
    K(ii-1,ii) = -k;
    K(ii,ii)   = 2*k;
    M(ii,ii)   = m;
end

for jj=1:n
    alpha_j=jj*pi/(n+1);
    omega2(jj,jj)=omega_0^2*2*(1-cos(alpha_j));
    for kk=1:n
       P(kk,jj)=2*sin(alpha_j*kk); 
    end
    norm_j=sqrt(P(:,jj)'*M*P(:,jj));
    P(:,jj)=P(:,jj)/norm_j;
end


figure
set(gca,'Fontsize',15)
hold on
plot(P(:,1)/max(P(:,1)),'r.-', 'Linewidth',2,'Markersize',15)
plot(P(:,2)/max(P(:,2)),'k.-', 'Linewidth',2,'Markersize',15)
plot(P(:,3)/max(P(:,3)),'m.-', 'Linewidth',2,'Markersize',15)
print -djpeg modes_equalized.jpg

figure
set(gca,'Fontsize',15)
hold on
plot(P(:,1),'r.-', 'Linewidth',2,'Markersize',15)
plot(P(:,2),'k.-', 'Linewidth',2,'Markersize',15)
plot(P(:,3),'m.-', 'Linewidth',2,'Markersize',15)
print -djpeg modes.jpg



F=zeros(n,1);
F(1)=1;
X_ref=K\F;

mu=P'*F;

for ii=1:n
   alpha_i=ii*pi/(n+1);    
   omega2(ii)=omega_0^2*2*(1-cos(alpha_i)); 
   q(ii,1)=mu(ii)/(omega2(ii,ii)); 
end


figure

set(gca,'Fontsize',15)

plot(X_ref)
hold on

nb_modes=modes_list(1);
plot(P(:,1:nb_modes)*q(1:nb_modes),'r.-', 'Linewidth',2,'Markersize',15)
nb_modes=modes_list(2);
plot(P(:,1:nb_modes)*q(1:nb_modes),'k.-', 'Linewidth',2,'Markersize',15)
nb_modes=modes_list(3);
plot(P(:,1:nb_modes)*q(1:nb_modes),'m.-', 'Linewidth',2,'Markersize',15)
nb_modes=modes_list(4);
plot(P(:,1:nb_modes)*q(1:nb_modes),'.-', 'Linewidth',2,'Markersize',15)

print -djpeg modal0.jpg

figure
set(gca,'Fontsize',15)
hold on
nb_modes=modes_list(1);
plot(P(:,1:nb_modes)*q(1:nb_modes)-X_ref,'r.-', 'Linewidth',2,'Markersize',15)
nb_modes=modes_list(2);
plot(P(:,1:nb_modes)*q(1:nb_modes)-X_ref,'k.-', 'Linewidth',2,'Markersize',15)
nb_modes=modes_list(3);
plot(P(:,1:nb_modes)*q(1:nb_modes)-X_ref,'m.-', 'Linewidth',2,'Markersize',15)
nb_modes=modes_list(4);
plot(P(:,1:nb_modes)*q(1:nb_modes)-X_ref,'.-', 'Linewidth',2,'Markersize',15)
print -djpeg errormodal0.jpg


X_modal=zeros(n,1);
for ii=1:n-1
    X_modal=X_modal+P(:,ii)*q(ii);
    error(ii)=(F'*F-F'*K*X_modal-X_modal'*K'*F+X_modal'*K^2*X_modal)/(F'*F);
end

