%%
%dados
clc
clear
a0 = 7.352460e+02;
a1 = 3.496580e-03;
a2 = -1.974370e-07;
a3 = 6.917050e-12;
a4 = -9.773650e-17;
b0 = 6.71633e+02;
b1 = 1.01738e-03;
b2 = -1.79972e-07;
b3 = 2.51328e-11;
b4 = 0.00000e+00;
k = 0.008633;  
beta=2.628;
qmin=200;
qmax=1692;
Xmin=5737.5; %5737.5 é múltiplo de 1% do reservatório, próximo ao 5773 previsto do exercício
Xmax=22950;
y = [740 606 500 412 433 509 721 1234 1770 1658 1480 1014];
d=1312*ones(1,12);
T=12;
X=[22950-(22950*0.01*27):22950*0.01:Xmax]; %estados do reservatório a passos de 1% do volume total
%X=[Xmin:22950*0.01:Xmax]; %estados do reservatório a passos de 1% do volume total

%%
%funções
custo = @(g) 0.02*g^2;
phi = @(x) a0+a1*x+a2*x^2+a3*x^3+a4*x^4;
theta = @(q) b0+b1*q+b2*q^2+b3*q^3+b4*q^4;
prod = @(x,q) k*q*(phi(x)-theta(q));
x_next = @(x,y,q) min(x+beta*(y-q),Xmax); 
%min(X_next,Xmax) devido qmax não ser um qvalid discreto para 
%casos em que ocorre vertimento
q_valid = @(x,x_next,y) y-(x_next-x)/beta; %obtem controle válidos para cada transição de estado
%%
c=inf*ones(size(X,2),T); %coluna T é o custo final de operação do reservatório
c(size(X,2),T)=0; %retira penalidade de custo para que o melhor custo final ocorra com reservatório cheio
q=qmin*ones(size(X,2),T);
for t=T:-1:1    
    for i=1:size(X,2) %i=estado no estágio atual (x)
        for j=1:size(X,2) %j=estado no estágio seguinte(x_next) 
            if(t==T) %ignorar último estágio
                continue
            end            
            q_loop=q_valid(X(i),X(j),y(t));   
            if(q_loop<qmin) 
                continue %controle infactível
            end
            if(q_loop>qmax) %excesso turbinado considero que seja vertimento
                q_loop=qmax;
            end
            g_loop=d(t)-prod(X(i),q_loop);
            if(g_loop<0) %excesso turbinado considero que seja vertimento
                g_loop=0;
            end            
            custo_loop=custo(g_loop);
            if(custo_loop+c(j,t+1)<c(i,t))
                c(i,t)=custo_loop+c(j,t+1);
                q(i,t)=q_loop;
            end    
        end
    end
end
%%
%plot pontos de cada estado possível
% for t=1:1:T
%     plot(t,X,'r.','MarkerSize',20);
%     hold on;
%     for i=1:size(X,2)
%         text(t,X(i)+100,strcat('c*=',num2str(fix(c(i,t)))));
%         text(t,X(i),strcat('q*=',num2str(fix(q(i,t)))));
%         if(t<T)
%             plot([t t+1], [X(i) x_next(X(i),y(t),q(i,t))]);
%         end;
%     end
% end
% axis([1-1 T+1 X(1)-1000 X(end)+1000]);
% xlabel('Estágio t');
% ylabel('Volume Reservatório [hm³]');

%%
%plot solução ótima
X_atual=22950;
for t=1:1:T
    plot(t,X_atual,'r.','MarkerSize',20);
    hold on;
    [row,col]=find(X==X_atual)
    text(t,X_atual+100,strcat('c*=',num2str(fix(c(col,t)))));
    text(t,X_atual,strcat('q*=',num2str(fix(q(col,t)))));
    if(t<T)
        plot([t t+1], [X_atual x_next(X_atual,y(t),q(col,t))], 'r');
        X_atual=x_next(X_atual,y(t),q(col,t))
    end;
end
axis([1-1 T+1 X(1)-1000 X(end)+1000]);
xlabel('Estágio t');
ylabel('Volume Reservatório [hm³]');
