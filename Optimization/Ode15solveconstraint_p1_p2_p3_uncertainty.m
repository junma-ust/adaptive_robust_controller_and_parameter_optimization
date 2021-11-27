close all
clear all
clc
global m1_ba m2_ba  r l g k1 k2  P  k3 kappa varepsilon 
m1_ba=2;
m2_ba=1;
r=0.1;
l=0.2;
g=9.8;

k1=1;
k2=0.1;
k3=0.1;
varepsilon=0.001;
P=1;

kappa_sanwei=[2.4458 1.4000 0.9603 0.5994 0.3759];

s0=[0 -0.6 0 -0.5 0.1]';%theta1 theta2 dtheta1 dtheta2 alpha 
ds0=[0 -0.5 0 0 0]';

for flag=1:1:5

kappa=kappa_sanwei(flag);

options=odeset; options.RelTol=1e-12;
[t,sout] = ode15i('Ode15sevroconstraint_p1_p2_p3_uncertainty',[0:1e-3:15],s0,ds0,options);
time=length(t);

for ka=1:time
    tte=t(ka);
%%%%%%%%%%%%%%%%%%adaptive robust control design%%%%%%%%%%%%%%%%%%%%%%
M_ba = [(m1_ba*3/2+m2_ba)*r^2 m2_ba*l*r*sin(sout(ka,2));m2_ba*l*r*sin(sout(ka,2)) 4*m2_ba*l^2/3];
Q1_ba=-m2_ba*l*r*sout(ka,4)^2*cos(sout(ka,2));
Q2_ba=-m2_ba*l*r*sout(ka,3)*sout(ka,4)*cos(sout(ka,2))+m2_ba*l*r*sout(ka,3)*sout(ka,4)^2*cos(sout(ka,2))-m2_ba*g*l*sout(ka,4)*cos(sout(ka,2));
Q_ba=[Q1_ba Q2_ba]';
B=[1 0]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%servo constraint%%%%%%%%%%%%%%%%%%%%%%%%
A=[0 1];
b=-sout(ka,3)-sout(ka,4)*sin(sout(ka,2))+0.1*cos(tte);
Adq_c=sout(ka,1)+sout(ka,4)+cos(sout(ka,2))-0.1*sin(tte);
Beta=Adq_c;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% control %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha_hat=sout(ka,5);
Pi_tilde=((sout(ka,3)^2+sout(ka,4)^2)^(1/2)+1)^2;
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

tau_w(ka)=p1(1)+p2(1)+p3(1);
error_0(ka)=Beta;
alpha_out(ka)=sout(ka,5);
end
color_bule = [0 0.4470 0.7410];
color_saffron_yellow = [0.8500, 0.3250, 0.0980];
color_yellow = [0.9290, 0.6940, 0.1250];
color_purple  = [0.4940, 0.1840, 0.5560];
color_green = [0.4660, 0.6740, 0.1880];
color_cyan = [0.3010, 0.7450, 0.9330];	
color_red = [0.6350, 0.0780, 0.1840];
color={color_bule,color_saffron_yellow,color_yellow,color_purple,color_green,color_cyan,color_red};

figure(11)
plot(t,tau_w,'LineWidth',2,'Color',color{flag});
hold on

figure(14)
plot(t,error_0,'LineWidth',2,'Color',color{flag});
hold on

figure(15)
plot(t,alpha_out,'LineWidth',2,'Color',color{flag});
hold on
end
figure(11)
grid on
legend({'$(\gamma_1,\gamma_2)=(100, 1)$','$(\gamma_1,\gamma_2)=(10, 1)$','$(\gamma_1,\gamma_2)=(1, 1)$','$(\gamma_1,\gamma_2)=(1, 10)$','$(\gamma_1,\gamma_2)=(1, 100)$'},'Interpreter','latex');
set(gca,'FontSize',32,'Fontname', 'Times New Roman');
xlabel('$t\ \rm{(s)}$','Interpreter','latex','FontSize',36,'Fontname', 'Times New Roman'); 
ylabel('$\tau \ \rm{(N \cdot m)}$','Interpreter','latex','FontSize',36,'Fontname', 'Times New Roman');

figure(14)
grid on
legend({'$(\gamma_1,\gamma_2)=(100, 1)$','$(\gamma_1,\gamma_2)=(10, 1)$','$(\gamma_1,\gamma_2)=(1, 1)$','$(\gamma_1,\gamma_2)=(1, 10)$','$(\gamma_1,\gamma_2)=(1, 100)$'},'Interpreter','latex');
set(gca,'FontSize',32,'Fontname', 'Times New Roman');
xlabel('$t\ \rm{(s)}$','Interpreter','latex','FontSize',36,'Fontname', 'Times New Roman'); 
ylabel('$\rm{Error} \ \theta_1+\dot{\theta}_2+cos(\theta_2)-0.1sin(t)$','Interpreter','latex','FontSize',36,'Fontname', 'Times New Roman');

figure(15)
grid on
legend({'$(\gamma_1,\gamma_2)=(100, 1)$','$(\gamma_1,\gamma_2)=(10, 1)$','$(\gamma_1,\gamma_2)=(1, 1)$','$(\gamma_1,\gamma_2)=(1, 10)$','$(\gamma_1,\gamma_2)=(1, 100)$'},'Interpreter','latex');
set(gca,'FontSize',32,'Fontname', 'Times New Roman');
xlabel('$t\ \rm{(s)}$','Interpreter','latex','FontSize',36,'Fontname', 'Times New Roman'); 
ylabel('$\hat{\delta}$','Interpreter','latex','FontSize',36,'Fontname', 'Times New Roman');
