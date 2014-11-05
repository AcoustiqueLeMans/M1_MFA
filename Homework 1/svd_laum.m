%   Program svd_laum.m
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
% This program performs a SVD of a picture of the LAUM. 
% It is associated to the course 
% "Mathematics for Acoustics and Acousticians"
%  given for students of the MSc Master Program in Acoustics and 
%  IMDEA students
% 


clear all 
close all
clc


disp('Importing LAUM Image')
image=double(imread('laum_color.jpg'));   
disp('SVD decomposition of the Matrix')

[U_R,S_R,V_R]=svd(image(:,:,1));
[U_G,S_G,V_G]=svd(image(:,:,2));
[U_B,S_B,V_B]=svd(image(:,:,3));

nb_images=rank(S_R);

i_temp=0;

for i_image=logspace(0,log10(nb_images),5)
  i_temp=i_temp+1;
  image_svd(:,:,1)=U_R(:,1:i_image)*S_R(1:i_image,1:i_image)*V_R(:,1:i_image)';
  image_svd(:,:,2)=U_G(:,1:i_image)*S_G(1:i_image,1:i_image)*V_G(:,1:i_image)';
  image_svd(:,:,3)=U_B(:,1:i_image)*S_B(1:i_image,1:i_image)*V_B(:,1:i_image)';
  imshow(uint8(image_svd))
  drawnow
  pause(1)
  print( ['svd_laum_ ' num2str(i_temp)],'-djpg')
end