% Algoritmo que genera las Curvas Teoricas de Q
% Parámetros del usuario para generar el diccionario

tamanoft=size(tfdist); % Calculo el tamaño del archivo de 
	X1=1; % X inicial
	X2=450; %X final
	Q1=1; % Q inicial
	Q2=150; %Q final es 150 porque se duplica
	
	X=X1:1:X2;
	%X=X*2*pi;
	Q=Q1:0.01:Q2;
	sizeX=size(X);
	sizeQ=size(Q);
Glim=20;

desvi=exp(-((0.23*Glim)+1.63));
	BetaQ=zeros(sizeX(2),sizeQ(2)); % Inicializa la matriz del Beta
	for i=1:sizeX(2) % Bucle de X
		for j=1:sizeQ(2) % Bucle de Q
		% Cálculo de los valores de teoricos de Q
			BetaQ(i,j)=exp(-(2*pi*i)/(2*Q(:,j)));
		end
	end

%Calculo de las variable usadas para determinar el control de amplitudes teoricos
	
	B1=BetaQ+desvi;
	B2=(BetaQ.^2)+desvi;
	Gan_Teorica=zeros(sizeX(2),sizeQ(2));
	Gan_Teorica=B1./B2;
	

%Figura de control de los datos de ganancia teoricos	
	
	numero=1;	
	figure;
	for i=1:sizeQ(2)
		%plot(Gan_Teorica(:,i));
		plot(X,Gan_Teorica(:,i));
		%axis([0,200,0,15]);
		title(['Ganancia Teorica ' num2str(numero)]);
		pause(0.001);
		numero=numero+0.01;
	end
	
	
savefile = 'Teorico_Q.mat'; % Crea un archivo para guardar el diccionario
save(savefile, 'BetaQ', 'Gan_Teorica') % Guarda el diccionario y otras variables en el archivo
clc
clear
	
	
	%%%%%Figura de control BetaQ
	figure
	numero=1;
	for i=1:sizeQ(2)
	plot(BetaQ(:,i));
	title(['Beta (X,Q) ' num2str(numero)]);
	pause(0.01);
	numero=numero+0.1;
	end
	
	
	%figure
	%for Q=1:Q2
	%plot(X,B1(:,Q),'b',X,B2(:,Q),'r');
	%plot(B2(:,Q));
	%axis([0,400,0,2])
	%title(['Beta (X,Q) ' num2str(Q)]);
	%pause(0.25);
	%end
	
	%for Q=1:Q2
	%plot(B2(:,Q));
	%axis([0,400,0,2])
	%title(['Beta (X,Q) ' num2str(Q)]);
	%pause(0.25);
	%end