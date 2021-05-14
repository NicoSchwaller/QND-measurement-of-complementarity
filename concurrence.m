function c=concurrence(rho)

try
Sigma = [0 0 0 -1 ; 0, 0, 1, 0 ; 0, 1, 0, 0 ; -1, 0, 0, 0];
V_p_rho=eig(rho);
R = rho*Sigma*rho'*Sigma;
V_p=eig(R);
V_p=sort(V_p);
arg1=0;
arg2=sqrt(abs(V_p(4)))-sqrt((abs(V_p(3))))-sqrt((abs(V_p(2))))-sqrt(abs(V_p(1)));
c=max(arg1,arg2);
catch
    c=NaN;
end


end