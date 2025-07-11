clc;clear all; close all;
Ts=0.01;KMAX=2500;T=Ts*KMAX;
m=.1;Fricc=0.1; long=0.6;g=9.8;M=.5;
TamanioFuente=12;
%pkg load signal; pkg load control; %Solo una vez
%Condiciones iniciales
alfa(1)=pi-.5; colorc='.g'; colorl='g';
alfa(1)=pi-.1; colorc='.r'; colorl='r';
% alfa(1)=pi-.8; colorc='.b';colorl='b';
%Versión linealizada en el equilibrio inestable. Sontag Pp 104.
%            estado=[p(i); p_p(i); alfa(i); omega(i)]
Mat_Ac=[0 1 0 0;0 -Fricc/M -m*g/M 0; 0 0 0 1; 0 Fricc/(long*M) g*(m+M)/(long*M) 0];
Mat_Bc=[0; 1/M; 0; -1/(long*M)];
% Equilibrio estable
Mat_A=[0 1 0 0;0 -Fricc/M -m*g/M 0; 0 0 0 1; 0 -Fricc/(long*M) -g*(m+M)/(long*M) 0];
Mat_B=[0; 1/M; 0; 1/(long*M)];
xOP=[0;0;pi;0];
Mat_C=[1 0 0 0;0 0 1 0]; %Mide delt y fi.
% Construcci?n del sistema ampliado
Mat_Aa=[Mat_A zeros(4,1);-Mat_C(1,:) 0];
Mat_Ba=[Mat_B;0];
Qa=diag([1e1 1e1 1e4 1e4 1e-3]);Ra=1e1 ;%1e4*eye(2); %q4=100000;q3=1000000;q2=10000; q1=1000;
Q=Qa(1:4,1:4);R=Ra;
Ha=[Mat_Aa -Mat_Ba*inv(Ra)*Mat_Ba'; -Qa -Mat_Aa'];
[n,va]=size(Ha);
[V,D]=eig(Ha);MX1X2=[];
for ii=1:n
    if real(D(ii,ii))<0
        MX1X2=[MX1X2 V(:,ii)];
    end
end
MX1=MX1X2(1:n/2,:); MX2=MX1X2(n/2+1:end,:);
Pa=real(MX2*inv(MX1));
Ka=inv(Ra)*Mat_Ba'*Pa;
eig(Mat_Aa-Mat_Ba*Ka)%  Polos del controlador
% break
K=Ka(1:4); KI=-Ka(5);
% abs(eig(Mat_A-Mat_B*K))
% break
%Observador:
%Repito con el Sistema Dual
Mat_Adual=Mat_A';
Mat_Bdual=Mat_C';
Mat_Cdual=Mat_B';
Qdual=diag([1e1 1 1e1 1]);Rdual= 1e-3*eye(2) ;%1e4*eye(2); %q4=100000;q3=1000000;q2=10000; q1=1000;
Ha=[Mat_Adual -Mat_Bdual*inv(Rdual)*Mat_Bdual'; -Qdual -Mat_Adual'];
[n,va]=size(Ha);
[V,D]=eig(Ha);MX1X2=[];
for ii=1:n
    if real(D(ii,ii))<0
        MX1X2=[MX1X2 V(:,ii)];
    end
end
MX1=MX1X2(1:n/2,:); MX2=MX1X2(n/2+1:end,:);
P=real(MX2*inv(MX1));
Ko=(inv(Rdual)*Mat_Bdual'*P)';
% break
t=0; x=[0;0;alfa(1);0];
% Position 0, ... , angle 0, ...
x=[0;0;0;0];
p(1)=x(1); p_p(1)=x(2); alfa(1)=x(3); omega(1)=x(4);
tl=(0:KMAX-1)*Ts;
x=[0;0;alfa(1);0];
x=[0;0;0;0];
h=Ts;ref=10;tita_pp=0;x_hat=[0;0;0;0];
um=0;%Zona muerta .1
for ki=1:KMAX
    estado=[p(ki); p_p(ki); alfa(ki); omega(ki)];
    Y_=Mat_C*estado;
    psi_p=ref-Mat_C(1,:)*estado;
    psi(ki+1)=psi(ki)+psi_p*h;    
    u(ki)=-K*(x_hat-xOP)+KI*psi(ki+1);
    % u(ki)=-K*(estado-xOP)+KI*psi(ki+1);
    acci(ki)=u(ki);
    %Zona muerta
    ui=u(ki);
    s2=u(ki);
    s22=0;
    if s2-um >= 0
        s22=max(0,s2-um);
    else
        if s2+um <= 0
            s22=min(0,s2+um);
        end
    end
    u(ki)=s22;
    %fin zona muerta
    p_pp=(1/(M+m))*(u(ki)-m*long*tita_pp*cos(alfa(ki))+m*long*omega(ki)^2*sin(alfa(ki))-Fricc*p_p(ki));
    tita_pp=(1/long)*(g*sin(alfa(ki))-p_pp*cos(alfa(ki)));
    p_p(ki+1)=p_p(ki)+h*p_pp;
    p(ki+1)=p(ki)+h*p_p(ki);
    omega(ki+1)=omega(ki)+h*tita_pp;
    alfa(ki+1)=alfa(ki)+h*omega(ki);
    %________OBSERVADOR__________%
    x_hatp=Mat_A*(x_hat-xOP)+Mat_B*u(ki)+Ko*(Y_-Mat_C*x_hat);
    x_hat=x_hat+h*x_hatp;
end
t=0:h:T;
figure(1);hold on;
subplot(3,2,1);plot(t,alfa,colorc);grid on;title('\phi_t','FontSize',TamanioFuente);hold on;
subplot(3,2,2);plot(t,omega,colorc);grid on;
title('$\dot{\phi_t}$','Interpreter','latex','FontSize',TamanioFuente);hold on;
subplot(3,2,3); plot(t,p,colorc);grid on;title('\delta_t','FontSize',TamanioFuente);hold on;
subplot(3,2,4);plot(t,p_p,colorc);grid on;title('$\dot{\delta_t}$','Interpreter','latex','FontSize',TamanioFuente);hold on;
subplot(3,1,3);plot(t(1:end-1),u,colorc,t(1:end-1),acci,'k');grid on;
title('Acción de control','FontSize',TamanioFuente);xlabel('Tiempo en Seg.','FontSize',TamanioFuente);hold on;

figure(2);hold on;
subplot(2,2,1);plot(alfa,omega,colorc);grid on;
xlabel('\phi_t','FontSize',TamanioFuente);
ylabel('$\dot{\phi_t}$','Interpreter','latex','Rotation',0,'FontSize',TamanioFuente);hold on;
subplot(2,2,2);plot(p,p_p,colorc);grid on;
xlabel('\delta_t','FontSize',TamanioFuente);hold on;
ylabel('$\dot{\delta_t}$','Interpreter','latex','Rotation',0,'FontSize',TamanioFuente);hold on;
% break
% figure(1);
% subplot(3,2,1);plot(tl,omegal,colorl,t,omega,colorc);grid on; title('Velocidad ángulo','FontSize',TamanioFuente);hold on;
% subplot(3,2,2);plot(tl,alfal,colorl,t,alfa,colorc);grid on;title('Ángulo','FontSize',TamanioFuente);hold on;
% % legend('Observador','Directo');legend('Boxoff')
% subplot(3,2,3);plot(tl,pl,colorl,t,p,colorc);grid on;title('Posición carro','FontSize',TamanioFuente);hold on;
% subplot(3,2,4);plot(tl,p_pl,colorl,t,p_p,colorc);grid on;title('Velocidad carro','FontSize',TamanioFuente);hold on;
% subplot(3,1,3);plot(tl,u_kl,colorl,t,u,colorc);grid on;title('Acción de control','FontSize',TamanioFuente);xlabel('Tiempo en Seg.','FontSize',TamanioFuente);hold on;
% figure(2);
% subplot(2,2,1);plot(alfal,omegal,colorl,alfa,omega,colorc);grid on;xlabel('Ángulo','FontSize',TamanioFuente);ylabel('Velocidad angular','FontSize',TamanioFuente);hold on;
% subplot(2,2,2);plot(pl,p_pl,colorl,p,p_p,colorc);grid on;xlabel('Posición carro','FontSize',TamanioFuente);ylabel('Velocidad carro','FontSize',TamanioFuente);hold on;
