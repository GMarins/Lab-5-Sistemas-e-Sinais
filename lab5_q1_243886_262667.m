%%%%%%%%%%%%%%%%%%LAB 5 - QUEST�O 1 %%%%%%%%%%%%%%%%
%Gabriel Marins da Costa - 243886
%Arthur Koucher - 262667
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Considere um sistema representado pela seguinte equa��o diferencial:
% y''' + 41 y'' + 360 y' + 900 y = 600 x' + 1200 x, onde y(0) = 2, y'(0) = 1
% e y''(0) = -0,05
% Baseando-se nisso, pede-se:

% (a) Calcular a solu��o desta equa��o diferencial, ou seja, calcular a 
% resposta natural, a resposta for�ada para um degrau unit�rio e a resposta 
% completa, como visto na �rea 1 da disciplina.
%
% (b) Determinar a transformada de Laplace do sistema (realizar os c�lculos
% a m�o ou usando a fun��o laplace).
%
% (c) Atrav�s da transformada de Laplace inversa (fun��o residue para
% fra��es parciais), obter as respostas natural, for�ada e completa do 
% sistema. Simular e comparar com as respostas obtidas no item (a).
%
% (d) Usando as fun��es de transfer�ncia, simular as respostas natural, 
% for�ada e completa usando os comandos step,impulse ou lsim. 
% Comparar com as respostas obtidas anteriormente.
%
% (e) Implementar um diagrama de blocos que representa a equa��o 
% diferencial dada no simulink (usando blocos somadores, integradores e 
% ganhos) e simule as respostas natural, for�ada e completa. Obtenha as 
% mesmas respostas usando blocos de fun��es de transfer�ncia.
% Conclus�o: Pudemos notar que o uso da Transformada de Laplace pode tornar
% a resolu��o de uma EDO muito mais simples e direta e tamb�m fornce um
% insight muito bom sobre o funcionamento do sistema, separado a resposta
% no dom�nio da frequ�ncia numa parte que depende de condi��es iniciais e
% em outra parte que depende apenas da entrada. As diversas simula��es sob
% circunst�ncias diferentes apresentaram as mesmas respostas, validando
% nossa resolu��o, o modelo utilizado e as ferramentas.


%Para o c�lculo anal�tico da equa��o, come�amos encontrando as ra�zes do
%polin�mio caracter�stico para a resposta natural
t = 0:0.001:3;
polos = roots([1 41 360 900]);

%A resposta natural tem a forma:
%yn(t) = Ae^p1 + Be^p2 + Ce^p3; onde pn s�o os polos
%As constantes s�o definidas resolvendo o sistema com as condi��es iniciais
%dadas
A = [1 1 1; -30 -6 -5; 900 36 25];
b = [2;1;-0.05];
n = A\b;
yn = n(1)*exp(polos(1)*t) + n(2)*exp(polos(2)*t) + n(3)*exp(polos(3)*t);
%Para o c�lculo da resposta for�ada, inspecionamos a entrada u(t) do
%sistema com a sua equa��o diferencial, considerando as condi��es iniciais
%nulas. Assim chegamos a uma resposta particular para a entrada u(t) dada e
%a uma resposta natural for�ada

%A solu��o particular � obtida por 
% 900*yp(t) = 1200
yp = 1200/900;

%Igualamos yf = yp+ynf �s condi��es iniciais
nf = A\[-4/3;0;600];
yf = yp + nf(1)*exp(polos(1)*t) + nf(2)*exp(polos(2)*t) + nf(3)*exp(polos(3)*t);

yc = yf + yn;

%Transformada feita a m�o utilizando a propriedade de diferencia��o -
%multiplica��o por s. Isolando Y(s) chegamos a dois membros do lado
%direito: um dependente das condi��es iniciais (resposta natural) e outro
%que n�o depende (resposta for�ada). 

fig1 = figure;
fig1.Name = 'Quest�o 1a)';

subplot(3,1,1);
plot(t,yn);
title('Resolu��o atrav�s do m�todo tradicional');
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

%Resposta natural obtida considerando x(t) = 0 e condi��es iniciais dadas 
%na propriedade de deriva��o da Transformada de Laplace

numerador = [2 83 760.95];
denominador = [1 41 360 900];
Hn = tf(numerador,denominador); %Fun��o de transfer�ncia dessa resposta
[z p k] = residue(numerador, denominador); %zn/(s-pn) + k 

syms s;
Yn = z./(s-p);
yn = ilaplace(Yn);

%Avaliamos yn[simb�lica] nos valores de t e convertemos para um resultado
%num�rico

yn = double(subs(yn)); %Subs avalia a fun��o e double converte para num�rico
yn = sum(yn,1);
%Resposta for�ada obtida considerando condi��es iniciais nulas

numerador = [600 1200];
denominador = [1 41 360 900 0];
[z p k] = residue(numerador, denominador);

syms s;
Yf = z./(s-p);
yf = ilaplace(Yf);

%Avaliamos yf[simb�lica] nos valores de t e convertemos para um resultado
%num�rico
yf = double(subs(yf));  %Subs avalia a fun��o e double converte para num�rico
yf = sum(yf,1);

yc = yn + yf; %Resposta completa

fig2 = figure;
fig2.Name = 'Quest�o 1c)';

subplot(3,1,1);
plot(t,yn);
title('Resolu��o atrav�s da Transformada de Laplace');
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
H = tf([600 1200],[1 41 360 900]); %Fun��o de Transfer�ncia do sistema
H_ss = ss(H); %Representa��o da TF no espa�o de estados

Ynat = impulse(Hn,t);
Yfor = step(H,t);
Ycomp = Ynat + Yfor;

fig3 = figure;
fig3.Name = 'Quest�o 1d)';

subplot(3,1,1);
plot(t,Ynat);
title('Resolu��o utilizando step e impulse');
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

%Simula��o com diagrama de blocos
%Para modelagem dos estados do sistema, isolamos o termo de maior grau e
%realizamos intetgra��es sucessivas. Pelo fato do sistema ser LTI, podemos
%reduzir o n�mero de blocos integradores at� o grau do sistema. Pelo mesmo
%fato, podemos manipular a ordem das opera��es mantendo a coes�o do
%sistema

sim('lab1_q1_blockdiagram');

fig4 = figure;
fig4.Name = 'Quest�o 1e) Diagrama de Blocos';

subplot(3,1,1);
plot(yn_blocos.Time,yn_blocos.Data);
title('Simula��o via Diagrama de Blocos');
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

%Para a simula��o utilizando a fun��o de transfer�ncia, nos aproveitamos do
%fato que Y(s) = H(s)X(s); H(s) �, portanto, autofun��o do sistema e,
%multiplicando a representa��o da entrada com a FT do sistema no dom�nio de
%Laplace, obtemos a transformada da sa�da

%De forma a excitar a entrada da resposta natural para atender �s
%especifica��es do programa, aplicamos um impulso instant�neo

fig5 = figure;
fig5.Name = 'Quest�o 1e) Simula��o com Fun��o de Transfer�ncia';

subplot(3,1,1);
plot(yn_transf.time,yn_transf.signals.values);
title('Simula��o via Fun��o de Transfer�ncia');
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