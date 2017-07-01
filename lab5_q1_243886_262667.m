%%%%%%%%%%%%%%%%%%LAB 5 - QUESTÃO 1 %%%%%%%%%%%%%%%%
%Gabriel Marins da Costa - 243886
%Arthur Koucher - 262667
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Considere um sistema representado pela seguinte equação diferencial:
% y''' + 41 y'' + 360 y' + 900 y = 600 x' + 1200 x, onde y(0) = 2, y'(0) = 1
% e y''(0) = -0,05
% Baseando-se nisso, pede-se:

% (a) Calcular a solução desta equação diferencial, ou seja, calcular a 
% resposta natural, a resposta forçada para um degrau unitário e a resposta 
% completa, como visto na área 1 da disciplina.
%
% (b) Determinar a transformada de Laplace do sistema (realizar os cálculos
% a mão ou usando a função laplace).
%
% (c) Através da transformada de Laplace inversa (função residue para
% frações parciais), obter as respostas natural, forçada e completa do 
% sistema. Simular e comparar com as respostas obtidas no item (a).
%
% (d) Usando as funções de transferência, simular as respostas natural, 
% forçada e completa usando os comandos step,impulse ou lsim. 
% Comparar com as respostas obtidas anteriormente.
%
% (e) Implementar um diagrama de blocos que representa a equação 
% diferencial dada no simulink (usando blocos somadores, integradores e 
% ganhos) e simule as respostas natural, forçada e completa. Obtenha as 
% mesmas respostas usando blocos de funções de transferência.
% Conclusão: Pudemos notar que o uso da Transformada de Laplace pode tornar
% a resolução de uma EDO muito mais simples e direta e também fornce um
% insight muito bom sobre o funcionamento do sistema, separado a resposta
% no domínio da frequência numa parte que depende de condições iniciais e
% em outra parte que depende apenas da entrada. As diversas simulações sob
% circunstâncias diferentes apresentaram as mesmas respostas, validando
% nossa resolução, o modelo utilizado e as ferramentas.


%Para o cálculo analítico da equação, começamos encontrando as raízes do
%polinômio característico para a resposta natural
t = 0:0.001:3;
polos = roots([1 41 360 900]);

%A resposta natural tem a forma:
%yn(t) = Ae^p1 + Be^p2 + Ce^p3; onde pn são os polos
%As constantes são definidas resolvendo o sistema com as condições iniciais
%dadas
A = [1 1 1; -30 -6 -5; 900 36 25];
b = [2;1;-0.05];
n = A\b;
yn = n(1)*exp(polos(1)*t) + n(2)*exp(polos(2)*t) + n(3)*exp(polos(3)*t);
%Para o cálculo da resposta forçada, inspecionamos a entrada u(t) do
%sistema com a sua equação diferencial, considerando as condições iniciais
%nulas. Assim chegamos a uma resposta particular para a entrada u(t) dada e
%a uma resposta natural forçada

%A solução particular é obtida por 
% 900*yp(t) = 1200
yp = 1200/900;

%Igualamos yf = yp+ynf às condições iniciais
nf = A\[-4/3;0;600];
yf = yp + nf(1)*exp(polos(1)*t) + nf(2)*exp(polos(2)*t) + nf(3)*exp(polos(3)*t);

yc = yf + yn;

%Transformada feita a mão utilizando a propriedade de diferenciação -
%multiplicação por s. Isolando Y(s) chegamos a dois membros do lado
%direito: um dependente das condições iniciais (resposta natural) e outro
%que não depende (resposta forçada). 

fig1 = figure;
fig1.Name = 'Questão 1a)';

subplot(3,1,1);
plot(t,yn);
title('Resolução através do método tradicional');
xlabel('t');
ylabel('${ y }_{ n }(t)$','Interpreter','LaTex');
axis([0 1.5 -0.5 3.5]);
grid on;

subplot(3,1,2);
plot(t,yf);
xlabel('t');
ylabel('${ y }_{ f }(t)$','Interpreter','LaTex');
axis([0 1.5 -0.5 3.5]);
grid on;

subplot(3,1,3);
plot(t,yc);
xlabel('t');
ylabel('${ y }_{ c }(t)$','Interpreter','LaTex');
axis([0 1.5 -0.5 3.5]);
grid on;

%%
%1c)
t = 0:0.001:3;

%Resposta natural obtida considerando x(t) = 0 e condições iniciais dadas 
%na propriedade de derivação da Transformada de Laplace

numerador = [2 83 760.95];
denominador = [1 41 360 900];
Hn = tf(numerador,denominador); %Função de transferência dessa resposta
[z p k] = residue(numerador, denominador); %zn/(s-pn) + k 

syms s;
Yn = z./(s-p);
yn = ilaplace(Yn);

%Avaliamos yn[simbólica] nos valores de t e convertemos para um resultado
%numérico

yn = double(subs(yn)); %Subs avalia a função e double converte para numérico
yn = sum(yn,1);
%Resposta forçada obtida considerando condições iniciais nulas

numerador = [600 1200];
denominador = [1 41 360 900 0];
[z p k] = residue(numerador, denominador);

syms s;
Yf = z./(s-p);
yf = ilaplace(Yf);

%Avaliamos yf[simbólica] nos valores de t e convertemos para um resultado
%numérico
yf = double(subs(yf));  %Subs avalia a função e double converte para numérico
yf = sum(yf,1);

yc = yn + yf; %Resposta completa

fig2 = figure;
fig2.Name = 'Questão 1c)';

subplot(3,1,1);
plot(t,yn);
title('Resolução através da Transformada de Laplace');
xlabel('t');
ylabel('${ y }_{ n }(t)$','Interpreter','LaTex');
axis([0 1.5 -0.5 3.5]);
grid on;

subplot(3,1,2);
plot(t,yf);
xlabel('t');
ylabel('${ y }_{ f }(t)$','Interpreter','LaTex');
axis([0 1.5 -0.5 3.5]);
grid on;

subplot(3,1,3);
plot(t,yc);
xlabel('t');
ylabel('${ y }_{ c }(t)$','Interpreter','LaTex');
axis([0 1.5 -0.5 3.5]);
grid on;

%%
%1d)
H = tf([600 1200],[1 41 360 900]); %Função de Transferência do sistema
H_ss = ss(H); %Representação da TF no espaço de estados

Ynat = impulse(Hn,t);
Yfor = step(H,t);
Ycomp = Ynat + Yfor;

fig3 = figure;
fig3.Name = 'Questão 1d)';

subplot(3,1,1);
plot(t,Ynat);
title('Resolução utilizando step e impulse');
xlabel('t');
ylabel('${ y }_{ n }(t)$','Interpreter','LaTex');
axis([0 1.5 -0.5 3.5]);
grid on;

subplot(3,1,2);
plot(t,Yfor);
xlabel('t');
ylabel('${ y }_{ f }(t)$','Interpreter','LaTex');
axis([0 1.5 -0.5 3.5]);
grid on;

subplot(3,1,3);
plot(t,Ycomp);
xlabel('t');
ylabel('${ y }_{ c }(t)$','Interpreter','LaTex');
axis([0 1.5 -0.5 3.5]);
grid on;

%%
%1e)

%Simulação com diagrama de blocos
%Para modelagem dos estados do sistema, isolamos o termo de maior grau e
%realizamos intetgrações sucessivas. Pelo fato do sistema ser LTI, podemos
%reduzir o número de blocos integradores até o grau do sistema. Pelo mesmo
%fato, podemos manipular a ordem das operações mantendo a coesão do
%sistema

sim('lab1_q1_blockdiagram');

fig4 = figure;
fig4.Name = 'Questão 1e) Diagrama de Blocos';

subplot(3,1,1);
plot(yn_blocos.Time,yn_blocos.Data);
title('Simulação via Diagrama de Blocos');
xlabel('t');
ylabel('${ y }_{ n }(t)$','Interpreter','LaTex');
axis([0 1.5 -0.5 3.5]);
grid on;

subplot(3,1,2);
plot(yf_blocos.Time,yf_blocos.Data);
xlabel('t');
ylabel('${ y }_{ f }(t)$','Interpreter','LaTex');
axis([0 1.5 -0.5 3.5]);
grid on;

subplot(3,1,3);
plot(yc_blocos.Time,yc_blocos.Data);
xlabel('t');
ylabel('${ y }_{ c }(t)$','Interpreter','LaTex');
axis([0 1.5 -0.5 3.5]);
grid on;

%Para a simulação utilizando a função de transferência, nos aproveitamos do
%fato que Y(s) = H(s)X(s); H(s) é, portanto, autofunção do sistema e,
%multiplicando a representação da entrada com a FT do sistema no domínio de
%Laplace, obtemos a transformada da saída

%De forma a excitar a entrada da resposta natural para atender às
%especificações do programa, aplicamos um impulso instantâneo

fig5 = figure;
fig5.Name = 'Questão 1e) Simulação com Função de Transferência';

subplot(3,1,1);
plot(yn_transf.time,yn_transf.signals.values);
title('Simulação via Função de Transferência');
xlabel('t');
ylabel('${ y }_{ n }(t)$','Interpreter','LaTex');
axis([0 1.5 -0.5 3.5]);
grid on;

subplot(3,1,2);
plot(yf_transf.time,yf_transf.signals.values);
xlabel('t');
ylabel('${ y }_{ f }(t)$','Interpreter','LaTex');
axis([0 1.5 -0.5 3.5]);
grid on;

subplot(3,1,3);
plot(yc_transf.time,yc_transf.signals.values);
xlabel('t');
ylabel('${ y }_{ c }(t)$','Interpreter','LaTex');
axis([0 1.5 -0.5 3.5]);
grid on;