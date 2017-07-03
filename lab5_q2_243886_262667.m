%%%%%%%%%%%%%%%%%%LAB 5 - QUESTÃO 2 %%%%%%%%%%%%%%%%
%Gabriel Marins da Costa - 243886
%Arthur Koucher - 262667
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Considere um sistema representado pela seguinte equação de diferenças
% (considerar o período de amostragem Ts = 0.1):
%           y[n] - 1.63y[n - 1] + 0.663y[n - 2] = x[n] - 0.8x[n - 2]
% onde y[-1] = 4 e demais condições iniciais são nulas. Baseado neste 
% sistema, pede-se:
% (a) Simular as respostas natural, forçada para um degrau unitário e 
% completa através da simulação da equação de recorrência.
%
% (b) Calcular a solução da equação de diferenças, ou seja, calcular a 
% resposta natural, a resposta forçada para um degrau unitário e a resposta
% completa, como visto na área 1 da disciplina.
% Comparar as respostas obtidas com as do item anterior.
% 
% (c) Determinar a transformada Z do sistema (realizar os cálculos à mão 
% ou através da função ztrans).
% 
% (d) Obter a transformada inversa do sistema para uma entrada x[n] = u[n] 
% e condições iniciais nulas através do método das divisões sucessivas 
% (também chamado "Série de Potências") e comparar com o resultado obtido 
% pela função step. Até qual coeficiente deve-se dividir o sistema para 
% termos uma boa aproximação da resposta?
% 
% (e) Obter a transformada inversa do sistema para entrada nula e condições
% iniciais diferentes de zero através do método das divisões sucessivas e 
% comparar com o resultado obtido pela função impulse. Até qual coeficiente
% deve-se dividir o sistema para termos uma boa aproximação da resposta?
%
% (f) Obter a transformada inversa do sistema para uma entrada x[n] = u[n] 
% através do método das frações parciais (comando residuez para frações 
% parciais). Calcule a resposta natural e a resposta forçada do sistema.
% Compare com o resultado do item anterior.
%
% (g) Implementar um diagrama de blocos que representa a equação de 
% diferenças dada no simulink (usando blocos somadores, de deslocamento e 
% ganhos) e simule as respostas natural, forçada e completa. Obtenha as 
% mesmas respostas usando blocos de funções de transferência.
% Conclusão: As diversas modelagens do mesmo sistema validaram todas as
% resoluções utilizadas. No caso específico de um sistema discreto, este é
% muito mais simples de se implementar computacionalmente devido à sua
% própria natureza recursiva. A resolução via Transformada Z se provou
% prática e rápida em comparação ao método analítico tradicional.

clear;
clc;

%2a
%A vantagem de se trabalhar com equações de diferenças é que sua modelagem
%por meio de programação é bastante direta e pode ser feita com um simples
%laço ou função de recorrência dadas as condições iniciais
%No MATLAB, podemos definir um vetor n de amostras e, com operações
%lógicas, ir preenchendo os vetores das respostas

n = -2:35;
%Para o cálculo da resposta natural, consideramos condições iniciais nulas
%y[n] = 1,63y[n-1] - 0,663y[n-2]; y[-1] = 4, y[-2] = 0
yn = 0 .* (n == -2); %Redundante, condição nula
yn = 4 .* (n == -1); %Aplicação da condição inicial utilizando expressão lógica
for i=0:length(n)
    yn(n == i) = 1.63*yn(n==i-1) - 0.663*yn(n==i-2);
end

%Resposta forçada considerando as entradas e condições iniciais nulas
% y[n]= 1.63*y[n-1] - 0.663*y[n-2] + u[n] - 0,8u[n-2]
yf = 0 .* (n == -2);
for i =0:length(n)
    yf(n == i) = 1.63*yf(n==i-1) - 0.663*yf(n==i-2) + 1.*(i>=0) - (8/10).*(i>=2);
end

yc = yn + yf;


fig1 = figure;
fig1.Name = 'Questão 2a)';

subplot(3,1,1);
stem(n,yn);
title('Resolução através da equação de recorrência');
xlabel('n');
ylabel('${ y }_{ n }[n]$','Interpreter','LaTex');
axis([0 35 0 10]);
grid on;

subplot(3,1,2);
stem(n,yf);
xlabel('n');
ylabel('${ y }_{ f }[n]$','Interpreter','LaTex');
axis([0 35 0 8]);
grid on;

subplot(3,1,3);
stem(n,yc);
xlabel('n');
ylabel('${ y }_{ c }[n]$','Interpreter','LaTex');
axis([0 35 0 15]);
grid on;

%%
%2b)
%Respostas analíticas calculadas a mão
n = -2:35;

yn = (41.2857*(0.85).^n - 34.7657*(0.78).^n)  .* (n>=0);
yf = (7.381*(0.85).^n - 12.4416*(0.78).^n + 6.0606) .*(n>=0);
yc = yn + yf;

H = tf([1 0 -0.8],[1 -1.63 0.663],1,'variable','z^-1'); %Função de transferência do sistema

fig2 = figure;
fig2.Name = 'Questão 2b)';

subplot(3,1,1);
stem(n,yn);
title('Resolução analítica');
xlabel('n');
ylabel('${ y }_{ n }[n]$','Interpreter','LaTex');
axis([0 35 0 10]);
grid on;

subplot(3,1,2);
stem(n,yf);
xlabel('n');
ylabel('${ y }_{ f }[n]$','Interpreter','LaTex');
axis([0 35 0 8]);
grid on;

subplot(3,1,3);
stem(n,yc);
xlabel('n');
ylabel('${ y }_{ c }[n]$','Interpreter','LaTex');
axis([0 35 0 15]);
grid on;

%%
%1d)
n = 0:35;

ystep = step(H,n);

%Para utilizar o método das divisões sucessivas, podemos usar a função
%deconv que realiza divisão polinomial

%Preparamos a representação de Laplace para a resposta forçada
num_forcada = [1 0 -8/10];
den_forcada = [1 -2.63 2.293 -0.663];

num_forcada = [num_forcada zeros(1,length(n))]; %Preparamos o numerador preenchendo-o com mais zeros
yforcada_coeficientes = deconv(num_forcada,den_forcada);


fig3 = figure;
fig3.Name = 'Questão 2d)';

subplot(2,1,1);
stem(n,yforcada_coeficientes);
title('Resposta forçada obtida através de divisões sucessivas');
xlabel('n');
ylabel('${ y }_{ f }[n]$','Interpreter','LaTex');
axis([0 35 0 8]);
grid on;

subplot(2,1,2);
title('Resposta forçada obtida através da excitação da função de transferência');
stem(n,ystep);
xlabel('n');
ylabel('${ y }_{ f }[n]$','Interpreter','LaTex');
axis([0 35 0 8]);
grid on;

%%
%1e)
%Realizamos o mesmo processo que na questão anterior, mas, dessa vez, com a
%transformada da resposta natural
n = 0:35;
num_natural = [6.52 -2.652];
den_natural = [1 -1.63 0.663];
Hn = tf(num_natural ,den_natural,1,'variable','z^-1');

num_natural = [num_natural zeros(1,length(n))];

ynatural_coeficientes = deconv(num_natural, den_natural);
yimpulse = impulse(Hn,n);

fig4 = figure;
fig4.Name = 'Questão 2d)';

subplot(2,1,1);
stem(n,ynatural_coeficientes);
title('Resposta natural obtida através de divisões sucessivas');
xlabel('n');
ylabel('${ y }_{ f }[n]$','Interpreter','LaTex');
axis([0 35 0 10]);
grid on;

subplot(2,1,2);
title('Resposta natural obtida através da excitação da função de transferência');
stem(n,yimpulse);
xlabel('n');
ylabel('${ y }_{ f }[n]$','Interpreter','LaTex');
axis([0 35 0 10]);
grid on;

%%
%1f)

num_natural = [6.52 -2.652];
den_natural = [1 -1.63 0.663];
num_forcada = [1 0 -8/10];
den_forcada = [1 -2.63 2.293 -0.663];

syms z;

[zn pn kn] = residuez(num_natural, den_natural);

Ynatural = zn./(z - pn) ;
ynatural = iztrans(Ynatural);

%Substituir os valores simbólicos pelos numéricos e somar os termos antes
%separados
ynatural = double(subs(ynatural));
ynatural = sum(ynatural,1);

[zf pf kf] = residuez(num_forcada, den_forcada);



Yforcada = zf./(z - pf);
yforcada = iztrans(Yforcada);

%Substituir os valores simbólicos pelos numéricos e somar os termos antes
%separados

yforcada = double(subs(yforcada));
yforcada = sum(yforcada,1);

ycompleta = ynatural + yforcada;

fig5 = figure;
fig5.Name = 'Questão 2f)';

subplot(3,1,1);
stem(n,ynatural);
title('Resolução utilizando frações parciais e Transformada Z inversa');
xlabel('n');
ylabel('${ y }_{ n }[n]$','Interpreter','LaTex');
axis([0 35 0 10]);
grid on;

subplot(3,1,2);
stem(n,yforcada);
xlabel('n');
ylabel('${ y }_{ f }[n]$','Interpreter','LaTex');
axis([0 35 0 8]);
grid on;

subplot(3,1,3);
stem(n,ycompleta);
xlabel('n');
ylabel('${ y }_{ c }[n]$','Interpreter','LaTex');
axis([0 35 0 15]);
grid on;

%%
%1g)
%A simulação via diagrama de blocos é bem direta para sistemas discretas.
%Sendo o sistema LTI, simplesmente isolamos y[n] e realizamos as operações
%de multiplicação por constante e atraso (integração em tempo discreta -
%blocos S).

%Podemos simular todas as três respostas usando apenas um diagrama; para
%isso definimos diferentes arrays para as diferentes condições de entrada
%ou condições iniciais
entrada = [0 1 1]; %0 Para R.Natural; 1 para R.Forçada e R.Completa
ys1 = [4 0 4]; %Condição inicial para y[-1]; 0 para R. Forçada
ys2 = [0 0 0]; %Condição inicial para y[-2]; 0 para todas

sim('lab5_q2_blockdiagram');

%A simulação retorna uma estrutura, onde yblocos.time é o vetor de n
%amostras e yblocos.signals.values é um array n x 3 em que cada coluna é
%uma das respostas do sistema

fig6 = figure;
fig6.Name = 'Questão 2g) Diagrama de Blocos';

subplot(3,1,1);
stem(yblocos.time,yblocos.signals.values(:,1));
title('Resolução utilizando diagrama de blocos');
xlabel('n');
ylabel('${ y }_{ n }[n]$','Interpreter','LaTex');
axis([0 35 0 10]);
grid on;

subplot(3,1,2);
stem(yblocos.time ,yblocos.signals.values(:,2));
xlabel('n');
ylabel('${ y }_{ f }[n]$','Interpreter','LaTex');
axis([0 35 0 8]);
grid on;

subplot(3,1,3);
stem(yblocos.time,yblocos.signals.values(:,3));
xlabel('n');
ylabel('${ y }_{ c }[n]$','Interpreter','LaTex');
axis([0 35 0 15]);
grid on;




% Para encontrarmos a resposta através de sua função de transferência,
% basta considerar o fato de que Y(z) = H(z)X(z) - ou seja, - H(z) é
% autofunção do sistema. O simulink possui um bloco específico para função
% de transferência, o que facilita a implementação


%
%Para simular a resposta utilizando blocos de função de transferência,
%basta utilizarmos os blocos disponíveis no simulink. Para a resposta
%natural, aplicamos um impulso na entrada do sistema a fim de excitá-lo
fig7 = figure;
fig7.Name = 'Questão 2g) Função de Transferência';

subplot(3,1,1);
stem(yftnatural.time,yftnatural.signals.values);
title('Resolução utilizando modelagem com função de transferência');
xlabel('n');
ylabel('${ y }_{ n }[n]$','Interpreter','LaTex');
axis([0 35 0 10]);
grid on;

subplot(3,1,2);
stem(yftforcada.time,yftforcada.signals.values);
xlabel('n');
ylabel('${ y }_{ f }[n]$','Interpreter','LaTex');
axis([0 35 0 8]);
grid on;

subplot(3,1,3);
stem(yftcompleta.time,yftcompleta.signals.values);
xlabel('n');
ylabel('${ y }_{ c }[n]$','Interpreter','LaTex');
axis([0 35 0 15]);
grid on;


