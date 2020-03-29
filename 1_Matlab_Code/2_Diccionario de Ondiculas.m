% ALGORITMO QUE GENERA EL DICCIONARIO DE OND�?CULAS MORLET
%Modificado para trabajar con una señal con muestreo de 2 ms y frecuencia maxima de 250 Hz
% Parámetros del usuario para generar el diccionario
	f1=2; % Frecuencia inicial
	f2=200; % Frecuencia final
	ph1=-180; % Fase inicial en grados
	ph2=180; % Fase final en grados
	ns=500; % Numero de muestras para una tasa de muestreo de 2 ms
	k=0.2; % Par�metro de ciclos de la ond�cula
	ts=0.5; % Extensi�n en tiempo de las ond�culas
	maxt=ns+1; % �ndice de la matriz
	f=f2-f1; % Diferencia de frecuencias
	ph=ph2-ph1; % Diferencia de fases
	rm=(ts*2)/(ns);
	t=(-ts:rm:ts); % Vector de tiempo de la ond�cula
	%T=(ts*2)/(ns); % Rata de muestreo
	frecuencia=1:f2;
	morletdic=zeros(f+1,ph+1,ns+1); % Inicializa la matriz del diccionario
    
		for m=1:f+1 % Bucle de frecuencias
			for n=1:ph+1 % Bucle de fases
				for o=1:ns+1 % Bucle de muestras en tiempo
					% C�lculo de los valores de la ond�cula de Morlet
					morletdic(m,n,o)=(exp((-1)*(t(o)^2)*((m-1)^2)*(log(2/k)))*exp(1i*2*pi*(m-1)*(t(o))))*exp(1i*(ph1+n-1)*pi/180);
				end
			end
		end

	savefile = 'diccionario_Qui_07_18_0_198.mat'; % Crea un archivo para guardar el diccionario
	save(savefile, 'morletdic','f1','f2','ph1','maxt','rm') % Guarda el diccionario y otras variables en el archivo
	clc
	clear
	

	figure(1)
	[x,y,z] = meshgrid(-180:1:180,0:1:200,-0.5:0.002:0.5); % Genera una grilla
	xslice = 180; yslice = [0,200]; zslice = 0.5; % Determina las caras a desplegar
	slice(x,y,z,real(morletdic),xslice,yslice,zslice) % Funci�n que realiza la gráfica
	colormap(gray) % Par�metro que establece la escala de color
	xlabel('Fase (grados)') % Etiqueta eje x
	ylabel('Frecuencia (hertz)') % Etiqueta eje y
	zlabel('Tiempo (segundos)') % Etiqueta eje z
	title('DICCIONARIO DE OND�?CULAS') % T�tulo de la gr�fica
	xlim([-180 180]) % limites eje fases
	ylim([0 200]) % limites eje frecuencias
	zlim([-1 1]) % limites eje tiempo