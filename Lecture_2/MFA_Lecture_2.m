%   Program MFA_L4.m
%
%   Copyright (C) 2015 LAUM UMR CNRS 6613 (France)
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

type_function=2

alpha=5.5*pi/2;
order=1+2*5;
is_Taylor=0;
is_Lagrange=1;
is_Legendre=1;

switch type_function
    case 1
        fonction = @(x) cos(alpha*x);
    case 2
        fonction = @(x) cos(alpha*x.^2);
end
 
x=linspace(-1,1,200);
f=fonction(x);

if type_function==2
    Gauss_points=compute_GaussLegendre_points(30);
    f_Gauss_points=fonction(Gauss_points.xi);
end

figure
set(gca,'Fontsize',25)
plot(x,f,'Linewidth',3,'Color','r')
hold on

if is_Taylor
    Taylor=ones(size(x));
    
    for kk=1:order
        coefficient=alpha^kk*real(1i^kk);
        Taylor=Taylor+coefficient*(x.^kk)/factorial(kk);
    end
    plot(x,Taylor,'Linewidth',3,'Color','k')
end

if is_Lagrange==1
x_Lagrange=linspace(-1,1,order);
f_Lagrange=fonction(x_Lagrange);
Lagrange=zeros(size(x));
for kk=1:order
    temp_Lagrange=ones(size(x));
    index=[1:kk-1 kk+1:order];
    for jj=index
        temp_Lagrange=temp_Lagrange.*(x-x_Lagrange(jj))/(x_Lagrange(kk)-x_Lagrange(jj));
    end
    temp_Lagrange=f_Lagrange(kk)*temp_Lagrange;
    Lagrange=Lagrange+temp_Lagrange;
end
plot(x,Lagrange,'Linewidth',3,'Color','b')
end

if is_Legendre
Legendre=zeros(size(x));
for kk=0:order-1
    L_k=Legendre_polynomial_normalized(kk);
    
    
    switch type_function
        case 1
            alpha_k=0;
            for jj=0:kk
                alpha_k=alpha_k+real(int_eialphax_xn(alpha,jj))*(L_k(jj+1));
            end
        case 2    
            % A numerical method to compute the integral
            alpha_k=sum(Gauss_points.w.*(f_Gauss_points.*evaluate_polynom_1D(L_k,Gauss_points.xi).'));
    end
    for jj=0:kk
        Legendre=Legendre+alpha_k*L_k(jj+1)*x.^(jj);
    end
end

plot(x,Legendre,'Linewidth',3,'Color','m')
end

ylim([-2 2])

figure
if is_Lagrange
plot(x,abs(f-Lagrange),'Linewidth',3,'Color','b')
hold on
end
if is_Legendre
plot(x,abs(f-Legendre),'Linewidth',3,'Color','m')
end
if is_Taylor
    plot(x,abs(f-Taylor),'Linewidth',3,'Color','k')
end
