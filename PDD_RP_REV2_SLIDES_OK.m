%%
%===============DEFINI��O DE DADOS======================
%=======================================================
clc
clear
%polinomio cota x volume
a0 = 7.352460e+02;
a1 = 3.496580e-03;
a2 = -1.974370e-07;
a3 = 6.917050e-12;
a4 = -9.773650e-17;
%polinomio de jusante
b0 = 6.71633e+02;
b1 = 1.01738e-03;
b2 = -1.79972e-07;
b3 = 2.51328e-11;
b4 = 0.00000e+00;
k = 0.008633;  
beta=2.6284; %convers�o m3/s -> hm3 (365/12)*24*60*60/1E6
qmin=200; %m�nima hist�rica deck NEWAVE � 102m3/s, mas a principio n�o ha limite
qmax=1692; %vaz�o dos dois cj de geradores 6*211+2*213 m3/s
Xmin=5737.5; %5737.5 � m�ltiplo de 1% do reservat�rio, pr�ximo ao 5733 dado deck NEWAVE
Xmax=22950;
Y_DATA = csvread('VAZOES_FURNAS_CSV.txt'); %HIST�RICO DECK 11/19
Y_DATA=Y_DATA(:,2:13);
Y_MLT=mean(Y_DATA(1:2019-1931,:));
% Y_MLT(1,:)=mean(Y_DATA(1:2019-1931,:));
% Y_MLT(2,:)=mean(Y_DATA(1:2009-1931,:));
% Y_MLT(3,:)=mean(Y_DATA(1:1999-1931,:));
% Y_MLT(4,:)=mean(Y_DATA(1:1989-1931,:));
% plot(1:12,Y_MLT);
% legend('MLT 1931-2018','MLT 1931-2008','MLT 1931-1998','MLT 1931-1988');
% title('Compara��o MLT');
% ylabel('Aflu�ncia [m�/s]');
% xlabel('M�s');
% ylim([0 2000])
Y_AUX=[];
for index=1:1:length(Y_DATA)
    Y_AUX=[Y_AUX Y_DATA(index,:)]; 
end

Y_HISTORICO=Y_AUX(1:(2019-1931)*12); %01/31 at� 12/18
%  
% y = [740 606 500 412 433 509 721 1234 1770 1658 1480 1014]; %MLT
y=[Y_MLT Y_MLT Y_MLT Y_MLT Y_MLT Y_MLT];
%%
%===============DEFINI��O DE FUN��ES========================
%===========================================================
custo = @(g) 0.02*g^2;
phi = @(x) a0+a1*x+a2*x^2+a3*x^3+a4*x^4;
theta = @(q) b0+b1*q+b2*q^2+b3*q^3+b4*q^4;
prod = @(x,q) k*q*(phi(x)-theta(q));
x_next = @(x,y,q) min(x+beta*(y-q),Xmax); 
defluencia_valida = @(x,x_next,y) y+(x/beta)-(x_next/beta); %obtem controle v�lidos para cada transi��o de estado
%%
%===============DEFINI��O DE VARI�VEIS======================
%===========================================================
T=length(y); 
X_otimo=NaN*zeros(1,T);
Q_otimo=NaN*zeros(1,T);
Vertimento_otimo=zeros(1,T);
G_otimo=zeros(1,T);
Custo_otimo=zeros(2,T);
d=1312*ones(1,T);
X=Xmin:Xmax*0.01:Xmax; %estados do reservat�rio a passos de 1% do volume total

c=inf*ones(size(X,2),T+1); %penaliza todos os �ltimos est�dos do est�gio T em inf
% c(size(X,2),T+1)=0; %retira penalidade de custo para que o melhor custo final ocorra com reservat�rio cheio
c(:,T+1)=0; %retira penalidade de custo para todos os estados do �ltimo est�gio
q=NaN*ones(size(X,2),T);
vertimento=zeros(size(X,2),T);
g=zeros(size(X,2),T);

%%
%===============ALGOR�TMO PDD BACKWARD======================
%===========================================================
tic
tstart=tic;
for t=T:-1:1    
    disp(strcat('t=',num2str(t)))
    for i=1:length(X) %i=estado no est�gio atual (x)
        for j=1:length(X) %j=estado no est�gio seguinte(x_next)   
            defluencia_loop=defluencia_valida(X(i),X(j),y(t));
            if(defluencia_loop>qmax) 
                 vert_loop=defluencia_loop-qmax;
                 q_loop=qmax;
            elseif(defluencia_loop<0) 
                continue %infactivel
            else
                vert_loop=0;
                q_loop=defluencia_loop;
            end
            g_loop=d(t)-prod(X(i),q_loop);
            while g_loop<0
                q_loop=q_loop-qmax*0.1;
                g_loop=d(t)-prod(X(i),q_loop);
                vert_loop=defluencia_loop-q_loop;
            end      
            custo_loop=custo(g_loop);
            if(custo_loop+c(j,t+1)<c(i,t))
                c(i,t)=custo_loop+c(j,t+1);
                q(i,t)=q_loop;
                vertimento(i,t)=vert_loop;
                g(i,t)=g_loop;
            end    
        end
    end
end
tempo_de_exec_PDD=toc(tstart)
%%
%===SIMULA��O A PARTIR DE X_atual e Y_DATA============
%===========================================================
X_atual=X(end); %partindo de reservat�rio cheio
% X_atual=X(1); %partindo de reservat�rio vazio
mes=1;
for t=1:1:length(Y_HISTORICO)   
    [val,idx]=min(abs(X-X_atual));
    X_otimo(t)=X_atual;
    Q_otimo(t)=q(idx,mes);
    Vertimento_otimo(t)=vertimento(idx,mes);
    G_otimo(t)=g(idx,mes);
    X_atual=x_next(X_atual,Y_HISTORICO(t),q(idx,mes));  
    mes=mes+1;
    if(mes==13)
        mes=1;
    end
end

% 
% X_atual=X(end); %partindo de reservat�rio vazio
% for t=1:1:length(y)   
%     [val,idx]=min(abs(X-X_atual));
%     X_otimo(1,t)=X_atual;
%     Q_otimo(1,t)=q(idx,t);
%     Vertimento_otimo(1,t)=vertimento(idx,t);
%     G_otimo(1,t)=g(idx,t);
%     X_atual=x_next(X_atual,y(t),q(idx,t));  
% end
% 
% y=[Y_MLT Y_MLT Y_MLT Y_MLT Y_MLT Y_MLT];
% c=inf*ones(size(X,2),T+1); %penaliza todos os �ltimos est�dos do est�gio T em inf
% % c(size(X,2),T+1)=0; %retira penalidade de custo para que o melhor custo final ocorra com reservat�rio cheio
% c(:,T+1)=0; %retira penalidade de custo para todos os estados do �ltimo est�gio
% q=NaN*ones(size(X,2),T);
% vertimento=zeros(size(X,2),T);
% g=zeros(size(X,2),T);
% T=length(y); 
% tic
% tstart=tic;
% for t=T:-1:1    
%     disp(strcat('t=',num2str(t)))
%     for i=1:length(X) %i=estado no est�gio atual (x)
%         for j=1:length(X) %j=estado no est�gio seguinte(x_next)   
%             defluencia_loop=defluencia_valida(X(i),X(j),y(t));
%             if(defluencia_loop>qmax) 
%                  vert_loop=defluencia_loop-qmax;
%                  q_loop=qmax;
%             elseif(defluencia_loop<0) 
%                 continue %infactivel
%             else
%                 vert_loop=0;
%                 q_loop=defluencia_loop;
%             end
%             g_loop=d(t)-prod(X(i),q_loop);
%             while g_loop<0
%                 q_loop=q_loop-qmax*0.1;
%                 g_loop=d(t)-prod(X(i),q_loop);
%                 vert_loop=defluencia_loop-q_loop;
%             end      
%             custo_loop=custo(g_loop);
%             if(custo_loop+c(j,t+1)<c(i,t))
%                 c(i,t)=custo_loop+c(j,t+1);
%                 q(i,t)=q_loop;
%                 vertimento(i,t)=vert_loop;
%                 g(i,t)=g_loop;
%             end    
%         end
%     end
% end
% tempo_de_exec_PDD=toc(tstart)
% 
% X_atual=X(end); %partindo de reservat�rio vazio
% for t=1:1:length(y)   
%     [val,idx]=min(abs(X-X_atual));
%     X_otimo(2,t)=X_atual;
%     Q_otimo(2,t)=q(idx,t);
%     Vertimento_otimo(2,t)=vertimento(idx,t);
%     G_otimo(2,t)=g(idx,t);
%     X_atual=x_next(X_atual,y(t),q(idx,t));  
% end


%%
%====LEITURA DE DADOS HYDROLAB PARA PLOT COMPARATIVO========
%===========================================================
X_DATA = csvread('VOL_PDD_HYDROLAB.txt');
X_otimo=[X_otimo;X_DATA.'];

Q_DATA = csvread('Q_PDD_HYDROLAB.txt');
Q_otimo=[Q_otimo;Q_DATA.'];

VERT_DATA = csvread('VERT_PDD_HYDROLAB.txt');
Vertimento_otimo=[Vertimento_otimo;VERT_DATA.'];

for m=1:length(X_otimo)
    PROD(m)=prod(X_otimo(1,m),Q_otimo(1,m));
    TOTALGER(m)=PROD(m)+G_otimo(m);
end

HPROD_DATA = csvread('HPROD_PDD_HYDROLAB.txt');
PROD=[PROD;HPROD_DATA.'];
GPROD_DATA = csvread('TPROD_PDD_HYDROLAB.txt');
G_otimo=[G_otimo;GPROD_DATA.'];
TOTALGER=[TOTALGER;PROD(2,:)+G_otimo(2,:)];

%C�LCULO CUSTO NA MESMA BASE PARA HYDROLAB
for p=length(X_otimo):-1:1
    if(p==length(X_otimo))
        Custo_otimo(1,p)=custo(G_otimo(1,p));
        Custo_otimo(2,p)=custo(G_otimo(1,p));
    else
        Custo_otimo(1,p)=Custo_otimo(1,p+1)+custo(G_otimo(1,p));
        Custo_otimo(2,p)=Custo_otimo(2,p+1)+custo(G_otimo(2,p));
    end
end

%===============SA�DA DE DADOS PARA GR�FICO=================
%===========================================================
T=length(Y_HISTORICO);
% x_axis=1:T;
x_axis=datetime(1930,12,1)+calmonths(1:T);

subplot(3,1,1);
[ax_volume nivel_reservatorio afluencia ] = plotyy(x_axis,X_otimo,x_axis,Y_HISTORICO);
set(nivel_reservatorio,'linestyle','-','marker','>');
set(ax_volume(1),'ylim',[5000 24000]);
ax_volume(1).YTick = [5000:(24000-5000)/5:24000];
title(strcat('Volume do reservat�rio e aflu�ncia por est�gio-',num2str(T),' meses'))
ylabel(ax_volume(1), 'Volume do Reservat�rio [hm�]');
xlabel('Est�gio t');% 
ylabel(ax_volume(2), 'Aflu�ncia [m�/s]');
set(afluencia,'linestyle','--');
legend('PDD','PDD Hydrolab','Aflu�ncia');
text(0+0.2,X_otimo(1,12)+250,strcat('Cota m�nima'));
text(0+0.2,X_otimo(1,1)+250,strcat('Cota m�xima'));
hold on
plot(x_axis,[ones(1,T)*Xmin;ones(1,T)*Xmax],'Color',[0 0 0],'LineStyle',':')
hold off

subplot(3,1,2);
[ax_turb turbinagem afluencia ] = plotyy(x_axis,Q_otimo,x_axis,Y_HISTORICO);
set(turbinagem,'linestyle','-','marker','>');
title(strcat('Turbinagem e aflu�ncia por est�gio-',num2str(T),' meses'));
ylabel(ax_turb(1), 'Turbinagem [m�/s]');
xlabel('Est�gio t');% 
ylabel(ax_turb(2), 'Aflu�ncia [m�/s]');
set(afluencia,'linestyle','--');
legend('PDD','PDD Hydrolab','Aflu�ncia');
hold on
plot(x_axis,[ones(1,T)*qmin;ones(1,T)*qmax],'Color',[0 0 0],'LineStyle',':')
hold off

subplot(3,1,3);
[ax_vert vert afluencia ] = plotyy(x_axis,Vertimento_otimo,x_axis,Y_HISTORICO);
set(vert,'linestyle','-','marker','>');
legend('PDD','PDD Hydrolab','Aflu�ncia');
title(strcat('Vertimento e aflu�ncia por est�gio-',num2str(T),' meses'));
ylabel(ax_vert(1), 'Vertimento [m�/s]');
xlabel('Est�gio t');% 
ylabel(ax_vert(2), 'Aflu�ncia [m�/s]');
set(afluencia,'linestyle','--');


linkaxes([ax_volume, ax_turb, ax_vert], 'x')


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% T=length(y);
% % x_axis=1:T;
% x_axis=1:T;
% 
% subplot(3,1,1);
% [ax_volume nivel_reservatorio afluencia ] = plotyy(x_axis,X_otimo,x_axis,y);
% set(nivel_reservatorio,'linestyle','-','marker','>');
% set(ax_volume(1),'ylim',[5000 24000]);
% ax_volume(1).YTick = [5000:(24000-5000)/5:24000];
% title('Volume do reservat�rio e aflu�ncia por est�gio')
% ylabel(ax_volume(1), 'Volume do Reservat�rio [hm�]');
% xlabel('Est�gio t');% 
% ylabel(ax_volume(2), 'Aflu�ncia [m�/s]');
% set(afluencia,'linestyle','--');
% legend('72 meses','60 meses','Aflu�ncia');
% hold on
% plot(x_axis,[ones(1,T)*Xmin;ones(1,T)*Xmax],'Color',[0 0 0],'LineStyle',':')
% hold off
% 
% subplot(3,1,2);
% [ax_turb turbinagem afluencia ] = plotyy(x_axis,Q_otimo,x_axis,y);
% set(turbinagem,'linestyle','-','marker','>');
% title('Turbinagem e aflu�ncia por est�gio');
% ylabel(ax_turb(1), 'Turbinagem [m�/s]');
% xlabel('Est�gio t');% 
% ylabel(ax_turb(2), 'Aflu�ncia [m�/s]');
% set(afluencia,'linestyle','--');
% legend('72 meses','60 meses','Aflu�ncia');
% hold on
% plot(x_axis,[ones(1,T)*qmin;ones(1,T)*qmax],'Color',[0 0 0],'LineStyle',':')
% hold off
% 
% subplot(3,1,3);
% [ax_vert vert afluencia ] = plotyy(x_axis,Vertimento_otimo,x_axis,y);
% set(vert,'linestyle','-','marker','>');
% legend('72 meses','60 meses','Aflu�ncia');
% title('Vertimento e aflu�ncia por est�gio');
% ylabel(ax_vert(1), 'Vertimento [m�/s]');
% xlabel('Est�gio t');% 
% ylabel(ax_vert(2), 'Aflu�ncia [m�/s]');
% set(afluencia,'linestyle','--');
% 
% 
% linkaxes([ax_volume, ax_turb, ax_vert], 'x')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure
subplot(3,1,1);
plot(x_axis,PROD);
legend('PDD Ger. Hidro.','PDD Hydrolab Ger. Hidro.');
title(strcat('Gera��o hidr�ulica por est�gio-',num2str(T),' meses'));
ylabel('Pot�ncia [MWmed]');
xlabel('Est�gio t'); 
ylim([0 1400])


subplot(3,1,2);
plot(x_axis,G_otimo);
legend('PDD Ger. Termo.','PDD Hydrolab Ger. Termo.');
title(strcat('Gera��o t�rmica por est�gio-',num2str(T),' meses'));
ylabel('Pot�ncia [MWmed]');
xlabel('Est�gio t');
ylim([0 1400])

subplot(3,1,3);
plot(x_axis,TOTALGER);
legend('PDD Ger. Total','PDD Hydrolab Ger. Total');
title(strcat('Gera��o total por est�gio-',num2str(T),' meses'));
ylabel('Pot�ncia [MWmed]');
xlabel('Est�gio t');
ylim([0 1400])

figure
plot(x_axis,Custo_otimo);
legend('PDD Custo Total','PDD Hydrolab Custo Total');
title(strcat('Custo total por est�gio-',num2str(T),' meses'));
ylabel('Custo [R$]');
xlabel('Est�gio t');

disp(strcat('Raz�o vertimento PDD/Hydrolab=',num2str(sum(Vertimento_otimo(1,:)/sum(Vertimento_otimo(2,:))))))
disp(strcat('Raz�o custo PDD/Hydrolab=',num2str(Custo_otimo(1,1)/Custo_otimo(2,1))))
disp(strcat('Raz�o produ��o hidroeletrica PDD/Hydrolab=',num2str(sum(PROD(1,:))/sum(PROD(2,:)))))


% sum(Vertimento_otimo(1,:))/sum(Vertimento_otimo(2,:))
% Custo_otimo(1,1)/Custo_otimo(2,1)
% sum(PROD(1,:))/sum(Q_otimo(1,:))
% sum(PROD(2,:))/sum(Q_otimo(2,:))
% 
% sum(PROD(1,:))/sum(PROD(2,:))
