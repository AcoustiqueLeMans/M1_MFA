%   Program int_eialphax_xn.m
%

function f=int_eialphax_xn(alpha,n)

% \int_{-1}^{1} e^{i \alpha x} x^n dx

if n==0
    f=(exp(1i*alpha)-exp(-1i*alpha))/(1i*alpha);
else
    f=(exp(1i*alpha)-(-1)^n*exp(-1i*alpha))/(1i*alpha);
    f=f-(n/(1i*alpha))*int_eialphax_xn(alpha,n-1);
end