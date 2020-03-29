%Data \Codigo-Tesis\QUI-07-17_GATHERS.sgy 
%CDP ESCOGIDOS
%Desde 2009 al 2200 (Comienzo, Pozo(RB-4)
%Desde 2700 hasta 2900
%Desde 3300 hasta 3600 Pozo(RB-53)

%IMPLEMENTACIÓN DEL ALGORITMO MATCHING PURSUIT DECOMPOSITION
% Previa la compilación de este algoritmo debe haber sido generado el diccionario de ondículas
% Se puede introducir una traza o una sección sísmica en forma de matriz
% Las filas deben ser las muestras en tiempo o profundidad
% Las columnas deben ser cada traza


clear all 							% Borra todas las variables del ambiente MATLAB
%whos -file diccionario_Qui_07_18.mat 	% Lee el archivo con el diccionario de ondículas
load diccionario_Qui_07_18.mat 			% Carga el archivo con el diccionario de ondículas
%load('F:\OneDrive\UNAL\Codigo_Efrain\Data_m\Quifa_18\Quifa_18_PSTM.mat') % Se carga el PSTM de la linea Quifa 18
load('F:\OneDrive\UNAL\Codigo_Efrain\Data_m\Quifa_18\Traza300.mat')
% Datos de entrada definidos por el usuario
%load QUI_07_17_corte_1.mat 	
%load QUI_07_17_corte_2.mat % Entrada de la sección sísmica como una matriz bidimensional
%load QUI_07_17_corte_3.mat % Entrada de la sección sísmica como una matriz bidimensional
%load QUI_07_17_corte_4.mat % Entrada de la sección sísmica como una matriz bidimensional
%load QUI_07_17_corte_5.mat % Entrada de la sección sísmica como una matriz bidimensional
%load QUI_07_17_corte_6.mat % Entrada de la sección sísmica como una matriz bidimensional
%load QUI_07_17_corte_7.mat % Entrada de la sección sísmica como una matriz bidimensional
%load QUI_07_17_corte_8.mat % Entrada de la sección sísmica como una matriz bidimensional
%load QUI_07_17_corte_9.mat % Entrada de la sección sísmica como una matriz bidimensional
%load QUI_07_17_corte_10.mat % Entrada de la sección sísmica como una matriz bidimensional

%Quifa corte 1 - 2-9-2015 - 2-12-2015
%Quifa corte 2 - 2-14-2015 - 2-16-2015
%Quifa corte 3 - 9/19/2015- - 9-23-2015 con 0.0001
%Quifa corte 3 - 9/19/2015- - 9-23-2015 con 0.0001
%Quifa corte 4 - 9/25/2015- - 9-28-2015 con 0.0001
%Quifa corte 5 - 10/04/2015- 10/6/2015  con 0.0001
%Quifa corte 6 - 10/07/2015 - 10-10-2015  con 0.0001
%Quifa corte 7 - 10-10-2015 - 10-13-2015 con 0.0001
%Quifa corte 8 - 10-14-2015 - con 0.0001
%load('QUI_07_17_corte_7.mat')
%load('diccionario_Qui_07_17.mat')

rm=0.001;
SeccionSismica=QUI_07_18;
SeccionSismica=traza;
%SeccionSismica=QUI_07_17_corte;
cutoff=0.001; 		% Valor de umbral de energía para que terminar el proceso 


MPDm=struct;%Crea estructura para almacenar los datos de la matriz MPD
FTdistribution=struct; %Crea estructura para almacenar los datos de la matriz FTdistribution
% Datos de entrada para el proceso en general. No es necesario modificar

SeccionSismica=SeccionSismica'; 	% Se transpone la secciÃ³n sÃ­smica
tamano=size(SeccionSismica); 		% Determina las dimensiones de la secciÃ³n sÃ­smica
num_trazas=tamano(1); 				% Se establece el nÃºmero de trazas
ts=0:rm:(tamano(2)-1)*rm; 			% Genera el vector de tiempo
mtst=zeros(tamano(1),tamano(2));
tptraza=zeros(tamano(1),1);
clc
barra = waitbar(0,'Initializing waitbar...');%Crea barra para llevar control grafico del proceso

	% Proceso repetitivo de descomposición de cada traza
	for m=1:num_trazas
		%MPD=zeros(10000,5); 												% Inicializa la matriz que almacena los valores de la descomposición
		tst=SeccionSismica(m,:); 											% Renombra cada traza a la vez como una nueva variable
		energia=energycalc(tst); 											% Invoca la función que calcula la energía inicial de la traza
		perc_energ=1;														% Se declara la constante de porcentaje de energía inicial = 100%
		cont=1; 															% Inicializa la variable de conteo de iteraciones
		bestondicula=1:maxt; 												% Se inicializa el vector que almacena la ondícula que mejor ajusta
		tst_min=min(tst); 													% Valor de amplitud mínimo en la traza de entrada actual
		tst_max=max(tst); 													% Valor de amplitud máximo en la traza de entrada actual
		mitadond=(maxt-1)/2; 												% Número de muestras de media ondícula
		modeledtrace=zeros(1,tamano(2));
		
		% Bucle de iteraciones hasta alcanzar el umbral de energía definido
		tic;
		while perc_energ>cutoff
			%Invoca la funcion que calcula los atributos complejos y almacena los parametros en MPD
			[MPD(cont,1),MPD(cont,2),MPD(cont,3),MPD(cont,4),MPD(cont,5)]=cmplxattrib(tst,rm,ts);
			% Extrae la ondícula del diccionario de ondículas
			% Suma la ondícula en cada iteración para obtener la traza modelada
			if MPD(cont,5) > 0
			
				for j=1:maxt
					bestondicula(j)=MPD(cont,3)*real(morletdic(MPD(cont,5)-(f1)+1,MPD(cont,4)-(ph1)+1,j));
				end
						
			else
			
				for j=1:maxt
					bestondicula(j)=MPD(cont,3)*real(morletdic(1,MPD(cont,4)-(ph1)+1,j));
				end
			
			end
		
			[tst,modeledtrace]=wavesubadd(tst,tamano(2),maxt,MPD(cont,2),bestondicula,modeledtrace); % Extrae la ondícula de la traza
			newenergy=energycalc(tst); % Calcula la energía de la nueva traza sin la ondícula actual
			perc_energ=newenergy/energia; % Porcentaje de energía de la nueva traza con respecto a la inicial
			cont=cont+1;%Incremento del contador
		end
		tptraza(m,1)=toc;

		mtst(m,:)=modeledtrace;
		
		MPDm.(['MPD_',num2str(m)])=MPD;
		FTdistribution.(['TFDist_',num2str(m)])=TimeFreqDistElmVSP(MPD(:,:), cont, ts, tamano(2),m);
		
		clearvars MPD;
		
		perc = (m*100)/num_trazas;
		waitbar(perc/100,barra,sprintf('%d%% along...',perc,m))

	end


	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
barra = waitbar(0,'Initializing waitbar...');%Crea barra para llevar control grafico del proceso
	for m=1:num_trazas
		MPD=MPDm.(['MPD_',num2str(m)]);
		MPDsize=size(MPD);
		cont=MPDsize(1);
		FTdistribution.(['TFDist_',num2str(m)])=TimeFreqDistElm(MPD(:,:), cont, ts, tamano(2),m);
		perc = (m*100)/num_trazas;
		waitbar(perc/100,barra,sprintf('%d%% along...',perc,m))

	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	for m=1:2402
		TFdist=FTdistribution.(['TFDist_',num2str(m)]);
		A(m,1)=sum(sum(TFdist));
	end
	
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	f=0:1:250;
	t=0.002:0.002:4;
	for m=1:2402
		TFdist=FTdistribution.(['TFDist_',num2str(m)]);
		imagesc(f,t,TFdist);					% Funcíon que realiza la gráfica como una imagen % Establece la escala de colores
		colormap(jet);
		ylim([0 2]) 
		ylabel('Tiempo (ms)');											% Etiqueta eje y
		xlabel('Frecuencia (hertz)');									% Etiqueta eje x
		zlabel('Amplitud (dB)');											% Etiqueta eje z
		title(['DISTRIBUCIÓN TIEMPO vs FRECUENCIA - TFdist ' num2str(m)]);	% Titulo de la gráfica xlim([0 200]) % Limites eje x
		pause(0.5);
	end
	
	
	
	
	
	