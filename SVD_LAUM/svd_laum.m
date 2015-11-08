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
image=double(imread('laum.jpg'));   
disp('SVD decomposition of the RGB matrices')

[U_R,S_R,V_R]=svd(image(:,:,1));
[U_G,S_G,V_G]=svd(image(:,:,2));
[U_B,S_B,V_B]=svd(image(:,:,3));

nb_images=10;


nb_sv=floor(logspace(0,log10(max([rank(S_R) rank(S_G) rank(S_B)])),nb_images));

image_svd=image;
image_svd(:,:,2)=0;
image_svd(:,:,3)=0;
imshow(uint8(image_svd))
print('laum_red','-djpeg')
image_svd=image;
image_svd(:,:,1)=0;
image_svd(:,:,3)=0;
imshow(uint8(image_svd))
print('laum_green','-djpeg')
image_svd=image;
image_svd(:,:,1)=0;
image_svd(:,:,2)=0;
imshow(uint8(image_svd))
print('laum_blue','-djpeg')

figure(1)
loglog(nonzeros(diag(S_R)),'r','Linewidth',5)
set(gca,'Fontsize',15)
hold on
loglog(nonzeros(diag(S_G)),'g','Linewidth',5)
loglog(nonzeros(diag(S_B)),'b','Linewidth',5)
legend('Red','Green','Blue')
xlabel('# of the singular value')
ylabel('Singular value')
print('laum_sv','-djpeg')

figure(3)

for ii=1:3

  image_svd(:,:,1)=U_R(:,ii)*S_R(ii,ii)*V_R(:,ii)';
  image_svd(:,:,2)=U_G(:,ii)*S_G(ii,ii)*V_G(:,ii)';
  image_svd(:,:,3)=U_B(:,ii)*S_B(ii,ii)*V_B(:,ii)';
  imshow(uint8(image_svd))
  title(['SVD Image of LAUM, Singular value #',num2str(ii)])
  drawnow
  pause(0.5)
  print( ['svd_laum_sv_ ' num2str(ii)],'-djpeg')
  
  error(ii)=sqrt(norm(image(:,:,1)-image_svd(:,:,1))^2+norm(image(:,:,2)-image_svd(:,:,2))^2+norm(image(:,:,3)-image_svd(:,:,3))^2);
  
end




for ii=1:nb_images

  image_svd(:,:,1)=U_R(:,1:nb_sv(ii))*S_R(1:nb_sv(ii),1:nb_sv(ii))*V_R(:,1:nb_sv(ii))';
  image_svd(:,:,2)=U_G(:,1:nb_sv(ii))*S_G(1:nb_sv(ii),1:nb_sv(ii))*V_G(:,1:nb_sv(ii))';
  image_svd(:,:,3)=U_B(:,1:nb_sv(ii))*S_B(1:nb_sv(ii),1:nb_sv(ii))*V_B(:,1:nb_sv(ii))';
  imshow(uint8(image_svd))
  title(['SVD Image of LAUM with ' num2str(nb_sv(ii)) ' Singular values'])
  drawnow
  pause(0.5)
  print( ['svd_laum_ ' num2str(ii)],'-djpeg')
  
  error(ii)=sqrt(norm(image(:,:,1)-image_svd(:,:,1))^2+norm(image(:,:,2)-image_svd(:,:,2))^2+norm(image(:,:,3)-image_svd(:,:,3))^2);
  error(ii)=error(ii)/sqrt(norm(image(:,:,1))^2+norm(image(:,:,2))^2+norm(image(:,:,3))^2);
end


figure
semilogy(nb_sv,error,'.-','Linewidth',5,'Markersize',15)
set(gca,'Fontsize',15)

xlabel('Number of Singular Values')
ylabel('Relative error')
print('laum_error','-djpeg')
