%Example extracted from: Sontag. Mathematical control theory 1998. Pag 104. http://www.sontaglab.org.
clc;clear all;close all;
m=.1;Fricc=0.1; long=0.6;g=9.8;M=.5;
h=0.001;tiempo=(120/h);tita_pp=0;
t=0:h:(tiempo)*h;
ref=100;
omega(1)=0; p_p(1)=0; u(1)=0; p(1)=0; i=1;indice=0;
%Versi?n linealizada en el equilibrio inestable. Sontag Pp 104.
% estado=[p(i); p_p(i); alfa(i); omega(i)]
Mat_A=[0 1 0 0;0 -Fricc/M -m*g/M 0; 0 0 0 1; 0 Fricc/(long*M) g*(m+M)/(long*M) 0];
Mat_B=[0; 1/M; 0; -1/(long*M)];
Mat_C=[1 0 0 0;0 0 1 0]; %La salida monovariable es posici?n
% Construcci?n del sistema ampliado
Mat_Aa=[Mat_A zeros(4,1);-Mat_C(1,:) 0];
Mat_Ba=[Mat_B;0];
% Q=np.diag([1e2, 1e2, 1e5, 1e5, 1]); R=np.diag([3e3]);
Qa=diag([1e2 1e2 1e5 1e5 1]);Ra=1e2 ;%1e4*eye(2); %q4=100000;q3=1000000;q2=10000; q1=1000;
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
% break
K=Ka(1:4); KI=-Ka(5);
eig(Mat_Aa-Mat_Ba*Ka)% break
%Observador:
%Repito con el Sistema Dual
Mat_Adual=Mat_A';
Mat_Bdual=Mat_C';
Mat_Cdual=Mat_B';
Qdual=diag([1 1 1 10000]);Rdual= 0.01*eye(2) ;%1e4*eye(2); %q4=100000;q3=1000000;q2=10000; q1=1000; 1 1 1 10000
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
x_hat=[0 0 0 0]'; Jn(1)=0;
T_CI=1;
Mat_datop=zeros(2,tiempo+1);
Mat_datop_p=zeros(2,tiempo+1);
Mat_datoalfa=zeros(2,tiempo+1);
Mat_datoomega=zeros(2,tiempo+1);
Mat_datou=zeros(2,tiempo+1);
hfig2 = figure(2); set(hfig2, 'Visible', 'on');hold on;
TamanioFuente=12;Ts=h;
um=0;
for ci =1:T_CI
    x_hat=[0 0 0 0]';Jn(1)=0;
    omega(1)=0; p_p(1)=0; u(1)=0; p(1)=0; i=1;indice=0;
    switch (ci)
        case 1
            alfa(1)=-pi;color='k';
        case 2
            alfa(1)=.6;color='r';
        case 3
            alfa(1)=-.6;color='g';
        case 4
            alfa(1)=pi;color='b';
    end
    colorc=[color '.'];
    Mat_datop(ci,1)=p(1);
    Mat_datop_p(ci,1)=p_p(1);
    Mat_datoalfa(ci,1)=alfa(1);
    Mat_datoomega(ci,1)=omega(1);
    Mat_datou(ci,1)=u(1);
    % eig(Mat_A-Mat_B*K)
    % break
    while(i<(tiempo+1))
        estado=[p(i); p_p(i); alfa(i); omega(i)];
        psi_p=ref-Mat_C(1,:)*estado;
        psi(i+1)=psi(i)+psi_p*h;
        %         u(i)=-K*estado+KI*psi(i+1);
        u11=-K*x_hat+KI*psi(i+1);
        ui=u11;
        s2=u11;
        s22=0;
        if s2-um >= 0
            s22=max(0,s2-um);
        else
            if s2+um <= 0
                s22=min(0,s2+um);
            end
        end
        u(i)=s22;
        
        % u(i)=max(-2000,u(i));u(i)=min(2000,u(i));
        p_pp=(1/(M+m))*(u(i)-m*long*tita_pp*cos(alfa(i))+m*long*omega(i)^2*sin(alfa(i))-Fricc*p_p(i));
        tita_pp=(1/long)*(g*sin(alfa(i))-p_pp*cos(alfa(i)));
        p_p(i+1)=p_p(i)+h*p_pp;
        p(i+1)=p(i)+h*p_p(i);
        omega(i+1)=omega(i)+h*tita_pp;
        alfa(i+1)=alfa(i)+h*omega(i);
        Y_=Mat_C*estado;
        y_sal(i)=Y_(1);
        Jn(i+1)=Jn(i)+estado'*Q(1:4,1:4)*estado+u(i)'*R*u(i);
        %________OBSERVADOR__________%
        y_sal_O=Mat_C*x_hat;
        x_hatp=Mat_A*x_hat+Mat_B*u(i)+Ko*(Y_-y_sal_O);
        x_hat=x_hat+h*x_hatp;
        i=i+1;
    end
    u(i)=-K*x_hat+KI*psi(i);
    Mat_datop(ci,:)=p;
    Mat_datop_p(ci,:)=p_p;
    Mat_datoalfa(ci,:)=alfa;
    Mat_datoomega(ci,:)=omega;
    Mat_datou(ci,:)=u;
    %     subplot(2,2,1);plot(alfa,omega,color);grid on;xlabel('Angle');ylabel('Angle velocity');hold on;
    %     subplot(2,2,2);plot(p,p_p,color);grid on;xlabel('Cart position');ylabel('Cart velocity');hold on;
    
    subplot(2,2,1);plot(alfa,omega,color);grid on;title('Ángulo','FontSize',TamanioFuente);
    xlabel('\phi_t','FontSize',TamanioFuente); 
    ylabel('$\dot{\phi_t}$','Interpreter','latex','Rotation',0,'FontSize',TamanioFuente);hold on;
    subplot(2,2,2);plot(p,p_p,color);grid on;xlabel('Posici?n carro','FontSize',TamanioFuente);hold on;
    
    
    ylabel('$\dot{\delta_t}$','Interpreter','latex','Rotation',0,'FontSize',TamanioFuente);hold on;
    xlabel('\delta_t','FontSize',TamanioFuente);
    title('Desplazamiento de carro','FontSize',TamanioFuente);
    subplot(2,1,2);semilogy(t,Jn,colorc);grid on;
    
    title('Modelo no lineal','FontSize',TamanioFuente);xlabel('Tiempo en Seg.','FontSize',TamanioFuente);ylabel('Funcional J','FontSize',TamanioFuente);hold on;
    ylabel('$J(x,u)$','Interpreter','latex','Rotation',0,'FontSize',TamanioFuente);hold on;
end
% print(hfig2,'Pendulo_fase','-dpng');
hfig1 = figure(1); set(hfig1, 'Visible', 'on');
subplot(3,2,1);plot(t,Mat_datoalfa);grid on;title('\phi_t','FontSize',TamanioFuente);hold on;
xlim([0, 60]);
subplot(3,2,2);plot(t,Mat_datoomega);grid on; title('$\dot{\phi_t}$','Interpreter','latex','FontSize',TamanioFuente);hold on; xlim([0, 60]);
subplot(3,2,3); plot(t,Mat_datop);grid on;title('\delta_t','FontSize',TamanioFuente);hold on;
xlim([0, 60]);
subplot(3,2,4);plot(t,Mat_datop_p);grid on;title('$\dot{\delta_t}$','Interpreter','latex','FontSize',TamanioFuente);hold on; xlim([0, 60]);
subplot(3,1,3);plot(t,Mat_datou);grid on;title('Acción de control','FontSize',TamanioFuente);xlabel('Tiempo en Seg.','FontSize',TamanioFuente);hold on;
xlim([0, 60]);
% print(hfig1,'Pendulo_t','-dpng');
% % save('Datos_Controlador.mat','-v7')
% break
% subplot(3,2,1);plot(t,Mat_datop);grid on; title('Cart position');hold on;
% subplot(3,2,2);plot(t,Mat_datoalfa);grid on;title('Angle');hold on;
% subplot(3,2,3); plot(t,Mat_datop_p);grid on;title('Cart velocity');hold on;
% subplot(3,2,4);plot(t,Mat_datoomega);grid on;title('Angle velocity');hold on;
% subplot(3,1,3);plot(t,Mat_datou);grid on;title('Force action');xlabel('Time Sec.');hold on;
% TamanioFuente=12;Ts=h;colorc=[color '.'];
% t=0:h:tiempo*(Ts);
% figure(1);hold on;
% subplot(3,2,1);plot(t,Mat_datoalfa,color);grid on;title('\phi_t','FontSize',TamanioFuente);hold on;
% subplot(3,2,2);plot(t,Mat_datoomega,color);grid on; title('$\dot{\phi_t}$','Interpreter','latex','FontSize',TamanioFuente);hold on;
% subplot(3,2,3); plot(t,Mat_datop,color);grid on;title('\delta_t','FontSize',TamanioFuente);hold on;
% subplot(3,2,4);plot(t,Mat_datop_p,color);grid on;title('$\dot{\delta_t}$','Interpreter','latex','FontSize',TamanioFuente);hold on;
% subplot(3,1,3);plot(t,u,color);grid on;title('Acci?n de control','FontSize',TamanioFuente);xlabel('Tiempo en Seg.','FontSize',TamanioFuente);hold on;
% figure(2);hold on;
% subplot(2,2,1);plot(Mat_datoalfa,Mat_datoomega,color);grid on;
% xlabel('\phi_t','FontSize',TamanioFuente);
% ylabel('$\dot{\phi_t}$','Interpreter','latex','Rotation',0,'FontSize',TamanioFuente);hold on;
% subplot(2,2,2);plot(Mat_datop,Mat_datop_p,color);grid on;xlabel('Posici?n carro','FontSize',TamanioFuente);hold on;
% ylabel('$\dot{\delta_t}$','Interpreter','latex','Rotation',0,'FontSize',TamanioFuente);hold on;
% xlabel('\delta_t','FontSize',TamanioFuente);
% subplot(2,1,2);semilogy(t,Jn,colorc);grid on;title('Modelo no lineal','FontSize',TamanioFuente);xlabel('Tiempo en Seg.','FontSize',TamanioFuente);ylabel('Funcional J','FontSize',TamanioFuente);hold on;
