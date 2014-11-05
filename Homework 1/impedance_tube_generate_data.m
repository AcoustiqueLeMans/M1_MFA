clear all
close all
clc

n=300;
x=linspace(-1,0,n);

f=440;
c=342;
k=2*pi*f/c;

noise=0.8
sigma=noise*rand(n,1);
error=sigma.*((rand(n,1)-0.5)+i*(rand(n,1)-0.5));

Inc=3
R=0.5*exp(j*pi/3)*Inc


P_theory=Inc*exp(-j*k*x')+R*Inc*exp(j*k*x');

maxP_theory=max(abs(P_theory));

P_theory_truncated=P_theory;

for ii=1:length(P_theory_truncated)
   if abs(P_theory_truncated(ii))<(0.6*maxP_theory)
      P_theory_truncated(ii)=0.6*maxP_theory*exp(j*2*pi*rand(1,1)); 
   end
end


P_measurement=P_theory_truncated+error;

fid_x=fopen('position.txt','w');
fid_theory=fopen('P_pure.txt','w');
fid_theory_truncated=fopen('P_background.txt','w');
fid_measurement=fopen('P_noisy.txt','w');

for ii=1:n
    fprintf(fid_x,'%1.15e\n',x(ii));
    fprintf(fid_theory,'%1.15e \t %1.15e\n',real(P_theory(ii)),imag(P_theory(ii)));
    fprintf(fid_theory_truncated,'%1.15e \t %1.15e\n',real(P_theory_truncated(ii)),imag(P_theory_truncated(ii)));
    fprintf(fid_measurement,'%1.15e \t %1.15e\n',real(P_measurement(ii)),imag(P_measurement(ii)));
end
fclose(fid_x);
fclose(fid_theory);
fclose(fid_theory_truncated);
fclose(fid_measurement);

figure(1)
hold on
plot(x,abs(P_theory))
plot(x,abs(P_measurement),'m')
figure(2)
hold on
plot(x,angle(P_theory))
plot(x,angle(P_measurement),'m')

