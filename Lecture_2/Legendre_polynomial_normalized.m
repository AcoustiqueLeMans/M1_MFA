function P_n=Legendre_polynomial_normalized(d)

switch d
    case 0
        P_n=[1]/sqrt(2);
    case 1
        P_n=[0 1]/sqrt(2/3);
    otherwise
        P_nm2=[1];
        P_nm1=[0 1];
        for n=2:d
            P_n=(((2*n)-1)/n)*multiply_polynom(P_nm1,[0 1])-((n-1)/n)*[P_nm2 0 0];
            P_nm2=P_nm1;
            P_nm1=P_n;
        end
        P_n=P_n/sqrt(2/(2*d+1));
end
