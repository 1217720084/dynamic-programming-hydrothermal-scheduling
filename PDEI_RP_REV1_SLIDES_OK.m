%%
%===============DEFINIÇÃO DE DADOS======================
%===========================================================
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
beta=2.6284; %conversão m3/s -> hm3 (365/12)*24*60*60/1E6
qmin=200; %mínima histórica deck NEWAVE é 102m3/s, mas a principio não ha limite
qmax=1692; %vazão dos dois cj de geradores 6*211+2*213 m3/s
Xmin=5737.5; %5737.5 é múltiplo de 1% do reservatório, próximo ao 5733 dado deck NEWAVE
Xmax=22950;
Y_DATA = csvread('VAZOES_FURNAS_CSV.txt'); %HISTÓRICO DECK 11/19 
Y_DATA=Y_DATA(:,2:13);
Y_AUX=[];
for index=1:1:length(Y_DATA)
    Y_AUX=[Y_AUX Y_DATA(index,:)]; 
end

% y=Y_AUX(1:(2019-1931)*12); %01/31 até 12/18
Y_HISTORICO=Y_AUX(1:(2019-1931)*12); %01/31 até 12/18
Y_MLT=mean(Y_DATA(1:2019-1931,:));

edges=0:500:4000;
for pos_mes=1:1:12
    [N,edges] = histcounts(Y_DATA(:,pos_mes), edges);
    prob(pos_mes,:)=N/length(Y_DATA(:,pos_mes));
end
y_edges=250:500:3750;
% histogram(Y_DATA(:,1),edges);
% title('Histograma afluência histórica Janeiro');
% xlabel('Blocos [m³/s]');
% ylabel('Frequência'); 
% plot(y_edges,prob(11,:))
% legend('jan','fev','mar','abr','mai','jun','jul','ago','set','out','nov','dez')
% [N,edges] = histcounts(Y_DATA(:,1));
% h = histogram(Y_DATA(:,1),edges);
% xgrid = linspace((edges(1)+edges(2))/2,(edges(end-1)+edges(end))/2, length(N));
% prob=N*/length(Y_DATA(:,1));
% line(xgrid,prob)
% trapz(xgrid,pdfEst)

% histfit(Y_DATA(:,1))

% %  
% y = [740 606 500 412 433 509 721 1234 1770 1658 1480 1014];
% y=[y y y y];
y=[Y_MLT Y_MLT Y_MLT Y_MLT Y_MLT Y_MLT Y_MLT Y_MLT Y_MLT Y_MLT];

%%
%===============DEFINIÇÃO DE FUNÇÕES========================
%===========================================================
custo = @(g) 0.02*g^2;
phi = @(x) a0+a1*x+a2*x^2+a3*x^3+a4*x^4;
theta = @(q) b0+b1*q+b2*q^2+b3*q^3+b4*q^4;
prod = @(x,q) k*q*(phi(x)-theta(q));
x_next = @(x,y,q) x+beta*(y-q); 
defluencia_valida = @(x,x_next,y) y+(x/beta)-(x_next/beta); %obtem controle válidos para cada transição de estado
%%
%===============DEFINIÇÃO DE VARIÁVEIS======================
%===========================================================
T=length(y); 
X_otimo=NaN*zeros(1,T);
Q_otimo=NaN*zeros(1,T);
Vertimento_otimo=zeros(1,T);
G_otimo=zeros(1,T);
Custo_otimo=zeros(2,T);
d=1312*ones(1,T);
X=Xmin:Xmax*0.01:Xmax; %estados do reservatório a passos de 1% do volume total
Q=qmin:qmax*0.01:qmax;
c=inf*ones(size(X,2),T+1); %penaliza todos os últimos estádos do estágio T em inf
% c(size(X,2),T+1)=0; %retira penalidade de custo para que o melhor custo final ocorra com reservatório cheio
c(:,T+1)=0; %retira penalidade de custo para todos os estados do último estágio
q=NaN*ones(size(X,2),T);
vertimento=zeros(size(X,2),T);
g=zeros(size(X,2),T);

%%
%===============ALGORÍTMO PDD BACKWARD======================
%===========================================================
tic
tstart=tic;
mes=12;
for t=T:-1:1    
    disp(strcat('t=',num2str(t)))    
    for i=1:length(X) %i=estado no estágio atual (x)
        for j=1:length(Q) %j=estado no estagio atual (q)
            flag_infactivel=0;
            g_loop=d(t)-prod(X(i),Q(j));
            if g_loop<0
                continue; %infactivel
            end
            custo_state=custo(g_loop);
            for k=1:1:length(y_edges)
                X_state=x_next(X(i),y_edges(k),Q(j));
                if  X_state < Xmin
                    flag_infactivel=1;
                    break; %infactivel
                end   
                [val,idx]=min(abs(X-X_state));
                custo_state=custo_state+c(idx,t+1)*prob(mes,k);
            end         
            if(custo_state<c(i,t) && flag_infactivel==0)
                c(i,t)=custo_state;
                q(i,t)=Q(j);
                g(i,t)=g_loop;
            end   
        end 
    end
    mes=mes-1
    if(mes==0)
        mes=12
    end
end

tempo_de_exec_PDEI=toc(tstart)
%%
%===CONSRTUÇÃO SOLUÇÃO ÓTIMA A PARTIR DE X_atual============
%===========================================================
X_atual=X(end); %partindo de reservatório cheio
mes=1;
for t=1:1:length(Y_HISTORICO)   
    [val,idx]=min(abs(X-X_atual));
    X_otimo(t)=X_atual;
    Q_otimo(t)=q(idx,mes);
    Vertimento_otimo(t)=vertimento(idx,mes);
    G_otimo(t)=g(idx,mes);
    X_atual=x_next(X_atual,Y_HISTORICO(t),q(idx,mes)); 
    if(X_atual>Xmax)
        Vertimento_otimo(t)=(X_atual-Xmax)/beta;
        X_atual=Xmax;
    end
    mes=mes+1;
    if(mes==13)
        mes=1;
    end
end

% 
% X_atual=X(end); %partindo de reservatório vazio
% for t=1:1:length(y)   
%     [val,idx]=min(abs(X-X_atual));
%     X_otimo(1,t)=X_atual;
%     Q_otimo(1,t)=q(idx,t);
%     Vertimento_otimo(1,t)=vertimento(idx,t);
%     G_otimo(1,t)=g(idx,t);
%     X_atual=x_next(X_atual,y(t),q(idx,t));  
%     if(X_atual>Xmax)
%         Vertimento_otimo(1,t)=(X_atual-Xmax)/beta;
%         X_atual=Xmax;
%     end
% end
% 
% 
% y=[Y_MLT Y_MLT Y_MLT Y_MLT Y_MLT ];
% T=length(y); 
% 
% c=inf*ones(size(X,2),T+1); %penaliza todos os últimos estádos do estágio T em inf
% % c(size(X,2),T+1)=0; %retira penalidade de custo para que o melhor custo final ocorra com reservatório cheio
% c(:,T+1)=0; %retira penalidade de custo para todos os estados do último estágio
% q=NaN*ones(size(X,2),T);
% vertimento=zeros(size(X,2),T);
% g=zeros(size(X,2),T);
% 
% tic
% tstart=tic;
% mes=12;
% for t=T:-1:1    
%     disp(strcat('t=',num2str(t)))    
%     for i=1:length(X) %i=estado no estágio atual (x)
%         for j=1:length(Q) %j=estado no estagio atual (q)
%             flag_infactivel=0;
%             g_loop=d(t)-prod(X(i),Q(j));
%             if g_loop<0
%                 continue; %infactivel
%             end
%             custo_state=custo(g_loop);
%             for k=1:1:length(y_edges)
%                 X_state=x_next(X(i),y_edges(k),Q(j));
%                 if  X_state < Xmin
%                     flag_infactivel=1;
%                     break; %infactivel
%                 end   
%                 [val,idx]=min(abs(X-X_state));
%                 custo_state=custo_state+c(idx,t+1)*prob(mes,k);
%             end         
%             if(custo_state<c(i,t) && flag_infactivel==0)
%                 c(i,t)=custo_state;
%                 q(i,t)=Q(j);
%                 g(i,t)=g_loop;
%             end   
%         end 
%     end
%     mes=mes-1
%     if(mes==0)
%         mes=12
%     end
% end
% 
% tempo_de_exec_PDEI=toc(tstart)
% 
% X_atual=X(end); %partindo de reservatório vazio
% for t=1:1:length(y)   
%     [val,idx]=min(abs(X-X_atual));
%     X_otimo(2,t)=X_atual;
%     Q_otimo(2,t)=q(idx,t);
%     Vertimento_otimo(2,t)=vertimento(idx,t);
%     G_otimo(2,t)=g(idx,t);
%     X_atual=x_next(X_atual,y(t),q(idx,t));  
%     if(X_atual>Xmax)
%         Vertimento_otimo(2,t)=(X_atual-Xmax)/beta;
%         X_atual=Xmax;
%     end
% end

%%
%====LEITURA DE DADOS HYDROLAB PARA PLOT COMPARATIVO========
%===========================================================
X_DATA = csvread('VOL_PDEI_HYDROLAB.txt');
X_otimo=[X_otimo;X_DATA.'];

Q_DATA = csvread('Q_PDEI_HYDROLAB.txt');
Q_otimo=[Q_otimo;Q_DATA.'];

VERT_DATA = csvread('VERT_PDEI_HYDROLAB.txt');
Vertimento_otimo=[Vertimento_otimo;VERT_DATA.'];

for m=1:length(X_otimo)
    PROD(m)=prod(X_otimo(1,m),Q_otimo(1,m));
    TOTALGER(m)=PROD(m)+G_otimo(m);
end

HPROD_DATA = csvread('HPROD_PDEI_HYDROLAB.txt');
PROD=[PROD;HPROD_DATA.'];
GPROD_DATA = csvread('TPROD_PDEI_HYDROLAB.txt');
G_otimo=[G_otimo;GPROD_DATA.'];
TOTALGER=[TOTALGER;PROD(2,:)+G_otimo(2,:)];

%CÁLCULO CUSTO NA MESMA BASE PARA HYDROLAB
for p=length(X_otimo):-1:1
    if(p==length(X_otimo))
        Custo_otimo(1,p)=custo(G_otimo(1,p));
        Custo_otimo(2,p)=custo(G_otimo(1,p));
    else
        Custo_otimo(1,p)=Custo_otimo(1,p+1)+custo(G_otimo(1,p));
        Custo_otimo(2,p)=Custo_otimo(2,p+1)+custo(G_otimo(2,p));
    end
end

%===============SAÍDA DE DADOS PARA GRÁFICO=================
%===========================================================
T=length(Y_HISTORICO);
% x_axis=1:T;
x_axis=datetime(1930,12,1)+calmonths(1:T);

subplot(3,1,1);
[ax_volume nivel_reservatorio afluencia ] = plotyy(x_axis,X_otimo,x_axis,Y_HISTORICO);
set(nivel_reservatorio,'linestyle','-','marker','>');
set(ax_volume(1),'ylim',[5000 24000]);
ax_volume(1).YTick = [5000:(24000-5000)/5:24000];
title(strcat('Volume do reservatório e afluência por estágio-',num2str(T),' meses'))
ylabel(ax_volume(1), 'Volume do Reservatório [hm³]');
xlabel('Estágio t');% 
ylabel(ax_volume(2), 'Afluência [m³/s]');
set(afluencia,'linestyle','--');
legend('PDEI','PDEI Hydrolab','Afluência');
text(0+0.2,X_otimo(1,12)+250,strcat('Cota mínima'));
text(0+0.2,X_otimo(1,1)+250,strcat('Cota máxima'));
hold on
plot(x_axis,[ones(1,T)*Xmin;ones(1,T)*Xmax],'Color',[0 0 0],'LineStyle',':')
hold off

subplot(3,1,2);
[ax_turb turbinagem afluencia ] = plotyy(x_axis,Q_otimo,x_axis,Y_HISTORICO);
set(turbinagem,'linestyle','-','marker','>');
title(strcat('Turbinagem e afluência por estágio-',num2str(T),' meses'));
ylabel(ax_turb(1), 'Turbinagem [m³/s]');
xlabel('Estágio t');% 
ylabel(ax_turb(2), 'Afluência [m³/s]');
set(afluencia,'linestyle','--');
legend('PDEI','PDEI Hydrolab','Afluência');
hold on
plot(x_axis,[ones(1,T)*qmin;ones(1,T)*qmax],'Color',[0 0 0],'LineStyle',':')
hold off

subplot(3,1,3);
[ax_vert vert afluencia ] = plotyy(x_axis,Vertimento_otimo,x_axis,Y_HISTORICO);
set(vert,'linestyle','-','marker','>');
legend('PDEI','PDEI Hydrolab','Afluência');
title(strcat('Vertimento e afluência por estágio-',num2str(T),' meses'));
ylabel(ax_vert(1), 'Vertimento [m³/s]');
xlabel('Estágio t');% 
ylabel(ax_vert(2), 'Afluência [m³/s]');
set(afluencia,'linestyle','--');


linkaxes([ax_volume, ax_turb, ax_vert], 'x')

% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% y=[Y_MLT Y_MLT Y_MLT Y_MLT Y_MLT Y_MLT Y_MLT Y_MLT Y_MLT Y_MLT];
% T=length(y);
% % x_axis=1:T;
% x_axis=1:T;
% 
% subplot(3,1,1);
% [ax_volume nivel_reservatorio afluencia ] = plotyy(x_axis,X_otimo,x_axis,y);
% set(nivel_reservatorio,'linestyle','-','marker','>');
% set(ax_volume(1),'ylim',[5000 24000]);
% ax_volume(1).YTick = [5000:(24000-5000)/5:24000];
% title('Volume do reservatório e afluência por estágio')
% ylabel(ax_volume(1), 'Volume do Reservatório [hm³]');
% xlabel('Estágio t');% 
% ylabel(ax_volume(2), 'Afluência [m³/s]');
% set(afluencia,'linestyle','--');
% legend('120 meses','60 meses','Afluência');
% hold on
% plot(x_axis,[ones(1,T)*Xmin;ones(1,T)*Xmax],'Color',[0 0 0],'LineStyle',':')
% hold off
% 
% subplot(3,1,2);
% [ax_turb turbinagem afluencia ] = plotyy(x_axis,Q_otimo,x_axis,y);
% set(turbinagem,'linestyle','-','marker','>');
% title('Turbinagem e afluência por estágio');
% ylabel(ax_turb(1), 'Turbinagem [m³/s]');
% xlabel('Estágio t');% 
% ylabel(ax_turb(2), 'Afluência [m³/s]');
% set(afluencia,'linestyle','--');
% legend('120 meses','60 meses','Afluência');
% hold on
% plot(x_axis,[ones(1,T)*qmin;ones(1,T)*qmax],'Color',[0 0 0],'LineStyle',':')
% hold off
% 
% subplot(3,1,3);
% [ax_vert vert afluencia ] = plotyy(x_axis,Vertimento_otimo,x_axis,y);
% set(vert,'linestyle','-','marker','>');
% legend('120 meses','60 meses','Afluência');
% title('Vertimento e afluência por estágio');
% ylabel(ax_vert(1), 'Vertimento [m³/s]');
% xlabel('Estágio t');% 
% ylabel(ax_vert(2), 'Afluência [m³/s]');
% set(afluencia,'linestyle','--');
% 
% 
% linkaxes([ax_volume, ax_turb, ax_vert], 'x')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



figure
subplot(3,1,1);
plot(x_axis,PROD);
legend('PDEI Ger. Hidro.','PDEI Hydrolab Ger. Hidro.');
title(strcat('Geração hidráulica por estágio-',num2str(T),' meses'));
ylabel('Potência [MWmed]');
xlabel('Estágio t'); 
ylim([0 1400])


subplot(3,1,2);
plot(x_axis,G_otimo);
legend('PDEI Ger. Termo.','PDEI Hydrolab Ger. Termo.');
title(strcat('Geração térmica por estágio-',num2str(T),' meses'));
ylabel('Potência [MWmed]');
xlabel('Estágio t');
ylim([0 1400])

subplot(3,1,3);
plot(x_axis,TOTALGER);
legend('PDEI Ger. Total','PDEI Hydrolab Ger. Total');
title(strcat('Geração total por estágio-',num2str(T),' meses'));
ylabel('Potência [MWmed]');
xlabel('Estágio t');
ylim([0 1400])

figure
plot(x_axis,Custo_otimo);
legend('PDEI Custo Total','PDEI Hydrolab Custo Total');
title(strcat('Custo total por estágio-',num2str(T),' meses'));
ylabel('Custo [R$]');
xlabel('Estágio t');

disp(strcat('Razão vertimento PDEI/Hydrolab=',num2str(sum(Vertimento_otimo(1,:)/sum(Vertimento_otimo(2,:))))))
disp(strcat('Razão custo PDEI/Hydrolab=',num2str(Custo_otimo(1,1)/Custo_otimo(2,1))))
disp(strcat('Razão produção hidroeletrica PDEI/Hydrolab=',num2str(sum(PROD(1,:))/sum(PROD(2,:)))))


% sum(Vertimento_otimo(1,:))/sum(Vertimento_otimo(2,:)))
% Custo_otimo(1,1)/Custo_otimo(2,1)
% sum(PROD(1,:))/sum(Q_otimo(1,:))
% sum(PROD(2,:))/sum(Q_otimo(2,:))
% 
% sum(PROD(1,:))/sum(PROD(2,:))

% csvwrite(strcat(num2str(T),'_VOL_PDEI.txt'),X_otimo(1,:).');
% csvwrite(strcat(num2str(T),'_TURB_PDEI.txt'),Q_otimo(1,:).');
% csvwrite(strcat(num2str(T),'_VERT_PDEI.txt'),Vertimento_otimo(1,:).');
% csvwrite(strcat(num2str(T),'_CUSTO_PDEI.txt'),Custo_otimo(1,:).');
% csvwrite(strcat(num2str(T),'_PRODT_PDEI.txt'),G_otimo(1,:).');
% csvwrite(strcat(num2str(T),'_TEMPO_EXEC_PDEI.txt'),tempo_de_exec_PDEI);