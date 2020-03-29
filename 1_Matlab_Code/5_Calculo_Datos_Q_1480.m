% ALGORITMO QUE CALCULA VALORES DE LA CURVA X=WT

tamanoM=size(tfdist);						%Identificar el tamaño del archivo de frecuencias
rm=0.002;									%Rata de muestreo para los datos es de 2 ms	
ts=0:rm:(tamanoM(1)-1)*rm; 					%Se genera un vector de tiempo
frecuencia=(0:0.5:200);						%Se genera un vector de frecuencias cada .5 Hz, siempre hasta 200 que e slo maximo del diccionario
tiempo=ts';									%Se traspone la matriz de timepo para crear vector Chi X
X=tiempo*frecuencia;						%Se crea vector Chi X
vector=unique(X);							%Se crea un vector de valores unicos de Chi
tamanoX=size(vector);						%Identifica cantidad de valores Chi
vector = vector(mod(vector, 0.002) == 0);	%Calculo valores multiplos a 0.002
vectorX=zeros(tamanoX(1),2);				%Crea un vectorX nuevo para almacenar las amplitudes por cada valor de Chi


%El siguiente ciclo for recorre l amatriz Chi y ubica cada posicion en la matriz de amplitud vs frecuencia
for i=1: tamanoX(1)								%Recorre la matriz de valores Chi
	Amplitud=0;									%Inicializa amplitud
	Temporal=0;									%Inicializa temporal
	cont=0;										%Inicializa contador
	value=vector(i);							%Inicializa valor del vector como primer dato CHi
	for j=1 : tamanoM(1)						%Recorre valores de la matriz de AvsF en direcion de las filas 
		for k=1 : tamanoM(2)					%Recorre valores de la matriz de AvsF en direcion de las columnas 
			if value == X(j,k)					%Si cumple con ser igual a l amatriz AvsF, guarda valor de ampltud y lo suma 
			Temporal=tfdist(j,k);
			Amplitud=Amplitud+(Temporal.^2); 
			cont=cont+1;
			end 
		end
	end
	vectorX(i,1)=Amplitud;
	vectorX(i,2)=cont;	
end




%Graficar del vector Chi
figure
plot(vectorX(:,1));							% Funcíon que realiza la gráfica como una imagen % Establece la escala de colores
ylabel('Amplitud')							% Etiqueta eje y
xlabel('Chi')								% Etiqueta eje x
title('Curva Chi') 							% Titulo de la gráfica xlim([0 200]) % Limites eje x



%Graficar del la Distribucion el campo de Isolineas
figure
imagesc(frecuencia,tiempo,tfdist);					% Funcíon que realiza la gráfica como una imagen % Establece la escala de colores
ylabel('Tiempo (ms)')						% Etiqueta eje y
xlabel('Frecuencia (hertz)')				% Etiqueta eje x
zlabel('Amplitud (dB)')						% Etiqueta eje z
title('DISTRIBUCIÓN TIEMPO vs FRECUENCIA') 	% Titulo de la gráfica xlim([0 200]) % Limites eje x
c = colorbar;


%Graficar el campo de Isolineas
figure
imagesc(frecuencia,tiempo,(X));
colormap(jet); % Establece la escala de colores
ylabel('Tiempo (ms)') % Etiqueta eje y
xlabel('Frecuencia Angular w (Radianes/Segundos)') % Etiqueta eje x
zlabel('Amplitud (dB)') % Etiqueta eje z
title('Distribucion Isolineas Tiempo VS Frecuencia') % Titulo de la gráfica
%xlim([0 frecuenciamax]) % Limites eje x
%ylim([0 tiempomax]) % Limites eje y

savefile = ['CurvaQdatos']; % Crea un archivo para guardar las curvas de los Xvalores %%%%num2str(nIso)
save(savefile, 'vectorX','vector') % Guarda el diccionario y otras variables en el archivo
clc
clear



















%Calculo de Q de los datos usando el Metodo de Compensacion de Amplitudes (WANG)

tamanox=size(Xwt_values);
[Xa,loc]=max(Xwt_values);

figure; 
plot(Xwt_values);      


XaNorm=(Xwt_values./Xa);

figure; 
plot(XaNorm);     

yXvalue=2*log(XaNorm);
plot(yXvalue); 


%Curva de ganancia controlada de los datos
Glim=20;
desvi=exp(-(0.23*Glim+1.63));
Gan_datos=zeros(tamanox(1),tamanox(2));

A=XaNorm+desvi;
B=(XaNorm.^2)+desvi;
Gan_datos=A./B;

for i=1:loc
Gan_datos(i,tamanox(2))=1;
end


figure; 
%plot(Gan_datos);   
plot(Gan_datos);  
%axis([0,200,0,12])
title('Ganancia Datos')


%RESTA DATOS CALCULADO vs TEORICO

a=Gan_datos;
b=Gan_Teorica;
c = zeros(size(b));
for i=1:size(b,2)
c(:,i) = b(:,i) - a;
end

valores=sum(c);
valores=abs(valores);

Qc=1:0.01:150;
figure; 
plot(Qc,valores);  
title('Ganancia Datos')

[Qcalculado,timeloc]=min(valores);






