clear all
close all
clc

Inc=3
R=0.5*exp(j*pi/3)*Inc

f=440;
c=342;
k=2*pi*f/c;

x=load('position.txt');
P_theory=load('P_pure.txt');
P_theory=P_theory(:,1)+i*P_theory(:,2);
P_theory_truncated=load('P_background.txt');
P_theory_truncated=P_theory_truncated(:,1)+i*P_theory_truncated(:,2);
P_measurement=load('P_noisy.txt');
P_measurement=P_measurement(:,1)+i*P_measurement(:,2);


figure(1)
hold on
plot(x,abs(P_theory))
plot(x,abs(P_theory_truncated),'k')
plot(x,abs(P_measurement),'m')
figure(2)
hold on
plot(x,angle(P_theory))
plot(x,angle(P_measurement),'m')
plot(x,angle(P_theory_truncated),'k')

M=[exp(-j*k*x) exp(j*k*x)];

X=(M'*M)\(M'*P_measurement)



P_fit=X(1)*exp(-j*k*x)+X(2)*exp(j*k*x);
figure(1)
hold on
plot(x,abs(P_fit),'r')
figure(2)
hold on
plot(x,angle(P_fit),'r')


for ii=1:100
    error_fit(ii)=norm([Inc;R]-X);
    W=diag(1./(abs(P_fit-P_measurement)));
    
    Mp=sqrt(W)*M;
    Pp=sqrt(W)*P_measurement;
    
    X=(Mp'*Mp)\(Mp'*Pp);
    
    P_fit=X(1)*exp(-j*k*x)+X(2)*exp(j*k*x);
    
    
end


figure(1)
hold on
plot(x,abs(P_fit),'c')
figure(2)
hold on
plot(x,angle(P_fit),'c')


figure
semilogy(error_fit,'.')
