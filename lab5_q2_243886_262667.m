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

