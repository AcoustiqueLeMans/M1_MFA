%   Program svd_laum_init.m
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


disp('Importing LAUM Image')
image=double(imread('laum.jpg'));   
disp('SVD decomposition of the Matrix')

[U_R,S_R,V_R]=svd(image(:,:,1));
[U_G,S_G,V_G]=svd(image(:,:,2));
[U_B,S_B,V_B]=svd(image(:,:,3));

for nb_sv=logspace(0,log10(917),30)

image_svd(:,:,1)=U_R(:,1:nb_sv)*S_R(1:nb_sv,1:nb_sv)*V_R(:,1:nb_sv)';
image_svd(:,:,2)=U_G(:,1:nb_sv)*S_G(1:nb_sv,1:nb_sv)*V_G(:,1:nb_sv)';
image_svd(:,:,3)=U_B(:,1:nb_sv)*S_B(1:nb_sv,1:nb_sv)*V_B(:,1:nb_sv)';
imshow(uint8(image_svd))

drawnow
pause(0.5)
end

