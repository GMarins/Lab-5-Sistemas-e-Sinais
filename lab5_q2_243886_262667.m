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

