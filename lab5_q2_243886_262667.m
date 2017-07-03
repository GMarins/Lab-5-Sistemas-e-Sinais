%%%%%%%%%%%%%%%%%%LAB 5 - QUEST�O 2 %%%%%%%%%%%%%%%%
%Gabriel Marins da Costa - 243886
%Arthur Koucher - 262667
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Considere um sistema representado pela seguinte equa��o de diferen�as
% (considerar o per�odo de amostragem Ts = 0.1):
%           y[n] - 1.63y[n - 1] + 0.663y[n - 2] = x[n] - 0.8x[n - 2]
% onde y[-1] = 4 e demais condi��es iniciais s�o nulas. Baseado neste 
% sistema, pede-se:
% (a) Simular as respostas natural, for�ada para um degrau unit�rio e 
% completa atrav�s da simula��o da equa��o de recorr�ncia.
%
% (b) Calcular a solu��o da equa��o de diferen�as, ou seja, calcular a 
% resposta natural, a resposta for�ada para um degrau unit�rio e a resposta
% completa, como visto na �rea 1 da disciplina.
% Comparar as respostas obtidas com as do item anterior.
% 
% (c) Determinar a transformada Z do sistema (realizar os c�lculos � m�o 
% ou atrav�s da fun��o ztrans).
% 
% (d) Obter a transformada inversa do sistema para uma entrada x[n] = u[n] 
% e condi��es iniciais nulas atrav�s do m�todo das divis�es sucessivas 
% (tamb�m chamado "S�rie de Pot�ncias") e comparar com o resultado obtido 
% pela fun��o step. At� qual coeficiente deve-se dividir o sistema para 
% termos uma boa aproxima��o da resposta?
% 
% (e) Obter a transformada inversa do sistema para entrada nula e condi��es
% iniciais diferentes de zero atrav�s do m�todo das divis�es sucessivas e 
% comparar com o resultado obtido pela fun��o impulse. At� qual coeficiente
% deve-se dividir o sistema para termos uma boa aproxima��o da resposta?
%
% (f) Obter a transformada inversa do sistema para uma entrada x[n] = u[n] 
% atrav�s do m�todo das fra��es parciais (comando residuez para fra��es 
% parciais). Calcule a resposta natural e a resposta for�ada do sistema.
% Compare com o resultado do item anterior.
%
% (g) Implementar um diagrama de blocos que representa a equa��o de 
% diferen�as dada no simulink (usando blocos somadores, de deslocamento e 
% ganhos) e simule as respostas natural, for�ada e completa. Obtenha as 
% mesmas respostas usando blocos de fun��es de transfer�ncia.
% Conclus�o: As diversas modelagens do mesmo sistema validaram todas as
% resolu��es utilizadas. No caso espec�fico de um sistema discreto, este �
% muito mais simples de se implementar computacionalmente devido � sua
% pr�pria natureza recursiva. A resolu��o via Transformada Z se provou
% pr�tica e r�pida em compara��o ao m�todo anal�tico tradicional.

clear;
clc;

%2a
%A vantagem de se trabalhar com equa��es de diferen�as � que sua modelagem
%por meio de programa��o � bastante direta e pode ser feita com um simples
%la�o ou fun��o de recorr�ncia dadas as condi��es iniciais
%No MATLAB, podemos definir um vetor n de amostras e, com opera��es
%l�gicas, ir preenchendo os vetores das respostas

n = -2:35;
%Para o c�lculo da resposta natural, consideramos condi��es iniciais nulas
%y[n] = 1,63y[n-1] - 0,663y[n-2]; y[-1] = 4, y[-2] = 0
yn = 0 .* (n == -2); %Redundante, condi��o nula
yn = 4 .* (n == -1); %Aplica��o da condi��o inicial utilizando express�o l�gica
for i=0:length(n)
    yn(n == i) = 1.63*yn(n==i-1) - 0.663*yn(n==i-2);
end

%Resposta for�ada considerando as entradas e condi��es iniciais nulas
% y[n]= 1.63*y[n-1] - 0.663*y[n-2] + u[n] - 0,8u[n-2]
yf = 0 .* (n == -2);
for i =0:length(n)
    yf(n == i) = 1.63*yf(n==i-1) - 0.663*yf(n==i-2) + 1.*(i>=0) - (8/10).*(i>=2);
end

yc = yn + yf;


fig1 = figure;
fig1.Name = 'Quest�o 2a)';

subplot(3,1,1);
stem(n,yn);
title('Resolu��o atrav�s da equa��o de recorr�ncia');
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
%Respostas anal�ticas calculadas a m�o
n = -2:35;

yn = (41.2857*(0.85).^n - 34.7657*(0.78).^n)  .* (n>=0);
yf = (7.381*(0.85).^n - 12.4416*(0.78).^n + 6.0606) .*(n>=0);
yc = yn + yf;

H = tf([1 0 -0.8],[1 -1.63 0.663],1,'variable','z^-1'); %Fun��o de transfer�ncia do sistema

fig2 = figure;
fig2.Name = 'Quest�o 2b)';

subplot(3,1,1);
stem(n,yn);
title('Resolu��o anal�tica');
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

%Para utilizar o m�todo das divis�es sucessivas, podemos usar a fun��o
%deconv que realiza divis�o polinomial

%Preparamos a representa��o de Laplace para a resposta for�ada
num_forcada = [1 0 -8/10];
den_forcada = [1 -2.63 2.293 -0.663];

num_forcada = [num_forcada zeros(1,length(n))]; %Preparamos o numerador preenchendo-o com mais zeros
yforcada_coeficientes = deconv(num_forcada,den_forcada);


fig3 = figure;
fig3.Name = 'Quest�o 2d)';

subplot(2,1,1);
stem(n,yforcada_coeficientes);
title('Resposta for�ada obtida atrav�s de divis�es sucessivas');
xlabel('n');
ylabel('${ y }_{ f }[n]$','Interpreter','LaTex');
axis([0 35 0 8]);
grid on;

subplot(2,1,2);
title('Resposta for�ada obtida atrav�s da excita��o da fun��o de transfer�ncia');
stem(n,ystep);
xlabel('n');
ylabel('${ y }_{ f }[n]$','Interpreter','LaTex');
axis([0 35 0 8]);
grid on;

%%
%1e)
%Realizamos o mesmo processo que na quest�o anterior, mas, dessa vez, com a
%transformada da resposta natural
n = 0:35;
num_natural = [6.52 -2.652];
den_natural = [1 -1.63 0.663];
Hn = tf(num_natural ,den_natural,1,'variable','z^-1');

num_natural = [num_natural zeros(1,length(n))];

ynatural_coeficientes = deconv(num_natural, den_natural);
yimpulse = impulse(Hn,n);

fig4 = figure;
fig4.Name = 'Quest�o 2d)';

subplot(2,1,1);
stem(n,ynatural_coeficientes);
title('Resposta natural obtida atrav�s de divis�es sucessivas');
xlabel('n');
ylabel('${ y }_{ f }[n]$','Interpreter','LaTex');
axis([0 35 0 10]);
grid on;

subplot(2,1,2);
title('Resposta natural obtida atrav�s da excita��o da fun��o de transfer�ncia');
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

%Substituir os valores simb�licos pelos num�ricos e somar os termos antes
%separados
ynatural = double(subs(ynatural));
ynatural = sum(ynatural,1);

[zf pf kf] = residuez(num_forcada, den_forcada);



Yforcada = zf./(z - pf);
yforcada = iztrans(Yforcada);

%Substituir os valores simb�licos pelos num�ricos e somar os termos antes
%separados

yforcada = double(subs(yforcada));
yforcada = sum(yforcada,1);

ycompleta = ynatural + yforcada;

fig5 = figure;
fig5.Name = 'Quest�o 2f)';

subplot(3,1,1);
stem(n,ynatural);
title('Resolu��o utilizando fra��es parciais e Transformada Z inversa');
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
%A simula��o via diagrama de blocos � bem direta para sistemas discretas.
%Sendo o sistema LTI, simplesmente isolamos y[n] e realizamos as opera��es
%de multiplica��o por constante e atraso (integra��o em tempo discreta -
%blocos S).

%Podemos simular todas as tr�s respostas usando apenas um diagrama; para
%isso definimos diferentes arrays para as diferentes condi��es de entrada
%ou condi��es iniciais
entrada = [0 1 1]; %0 Para R.Natural; 1 para R.For�ada e R.Completa
ys1 = [4 0 4]; %Condi��o inicial para y[-1]; 0 para R. For�ada
ys2 = [0 0 0]; %Condi��o inicial para y[-2]; 0 para todas

sim('lab5_q2_blockdiagram');

%A simula��o retorna uma estrutura, onde yblocos.time � o vetor de n
%amostras e yblocos.signals.values � um array n x 3 em que cada coluna �
%uma das respostas do sistema

fig6 = figure;
fig6.Name = 'Quest�o 2g) Diagrama de Blocos';

subplot(3,1,1);
stem(yblocos.time,yblocos.signals.values(:,1));
title('Resolu��o utilizando diagrama de blocos');
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




% Para encontrarmos a resposta atrav�s de sua fun��o de transfer�ncia,
% basta considerar o fato de que Y(z) = H(z)X(z) - ou seja, - H(z) �
% autofun��o do sistema. O simulink possui um bloco espec�fico para fun��o
% de transfer�ncia, o que facilita a implementa��o


%
%Para simular a resposta utilizando blocos de fun��o de transfer�ncia,
%basta utilizarmos os blocos dispon�veis no simulink. Para a resposta
%natural, aplicamos um impulso na entrada do sistema a fim de excit�-lo
fig7 = figure;
fig7.Name = 'Quest�o 2g) Fun��o de Transfer�ncia';

subplot(3,1,1);
stem(yftnatural.time,yftnatural.signals.values);
title('Resolu��o utilizando modelagem com fun��o de transfer�ncia');
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


