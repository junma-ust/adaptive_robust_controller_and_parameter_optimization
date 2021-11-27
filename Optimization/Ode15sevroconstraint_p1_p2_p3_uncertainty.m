function f=Ode15sevroconstraint_p1_p2_p3_uncertainty(t,s,ds)
global m1_ba m2_ba  r l g k1 k2 P  k3 kappa varepsilon 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%model%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M_ba = [(m1_ba*3/2+m2_ba)*r^2 m2_ba*l*r*sin(s(2));m2_ba*l*r*sin(s(2)) 4*m2_ba*l^2/3];
Q1_ba=-m2_ba*l*r*s(4)^2*cos(s(2));
Q2_ba=-m2_ba*l*r*s(3)*s(4)*cos(s(2))+m2_ba*l*r*s(3)*s(4)^2*cos(s(2))-m2_ba*g*l*s(4)*cos(s(2));
Q_ba=[Q1_ba Q2_ba]';
B=[1 0]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% model with uncertainty %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m1=m1_ba+0.2*m1_ba*cos(t);
m2=m2_ba+0.4*m2_ba*sin(t);

M = [(m1*3/2+m2)*r^2 m2*l*r*sin(s(2));m2*l*r*sin(s(2)) 4*m2*l^2/3];
Q1=-m2*l*r*s(4)^2*cos(s(2));
Q2=-m2*l*r*s(3)*s(4)*cos(s(2))+m2*l*r*s(3)*s(4)^2*cos(s(2))-m2*g*l*s(4)*cos(s(2));
Q=[Q1 Q2]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%servo constraint%%%%%%%%%%%%%%%%%%%%%%%%%%%

A=[0 1];
b=-s(3)-s(4)*sin(s(2))+0.1*cos(t);
Adq_c=s(1)+s(4)+cos(s(2))-0.1*sin(t);
Beta=Adq_c;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% control %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alpha_hat=s(5);
Pi_tilde=((s(3)^2+s(4)^2)^(1/2)+1)^2;
dalpha_hat=kappa*(k1*Pi_tilde*norm(Beta)-(k2+k3*exp(-norm(Beta)))*alpha_hat);
Pi=alpha_hat*Pi_tilde;

mu=Beta*Pi;

if norm(mu)>varepsilon
    gamma=1/norm(mu);
else
    gamma=1/varepsilon;
end

p1=pinv(A*M_ba^(-1)*B)*(b-A*M_ba^-1*Q_ba);
p2=-kappa*(A*M_ba^(-1)*B)'*((A*M_ba^(-1)*B)*(A*M_ba^(-1)*B)')^(-1)*P^(-1)*Beta;
p3=-(A*M_ba^(-1)*B)'*((A*M_ba^(-1)*B)*(A*M_ba^(-1)*B)')^(-1)*P^(-1)*gamma*mu*Pi;

ddq=M^(-1)*Q+M^(-1)*B*(p1+p2+p3);

f=[ds(1)-s(3);
   ds(2)-s(4);
   ds(3)-ddq(1);
   ds(4)-ddq(2);
   ds(5)-dalpha_hat
   ]
