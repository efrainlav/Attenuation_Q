%Leer los segy de las lineas QUI17 - QUI18 y cargarlo en Matlab
%de tiempo hasta 4000 ms

%Se realiza el llamado del segy
QUI_07_17_100=read_segy_file('J:\Estudio\Maestria\Tesis\DATA\QUI-07-17_GATHERS.sgy',{'traces',1:100});%Selecciona de l atraza 1 a la 100 y crea un archivo segy en formato matlab
QUI_07_17=read_segy_file('J:\Estudio\Maestria\Tesis\DATA\QUI-07-17_GATHERS.sgy');%carga el seg completo con los header a matlab
QUI_07_18=read_segy_file('J:\Estudio\Maestria\Tesis\DATA\QUI-07-18_GATHERS.sgy');%carga el seg completo con los header a matlab
QUI_07_17_Traza1000=read_segy_file('J:\Estudio\Maestria\Tesis\DATA\QUI-07-17_GATHERS.sgy',{'traces',1000});%Selecciona solo la traza 1000

%Ploteamos para verificar la señal
s_wplot(QUI_07_17_Traza1000);

%Se crea un segy de solo la traza 1000 
write_segy_file(QUI_07_17_Traza1000); % Se llamo Traza_1000


%Se realiza el llamado de nuevo del segy pero creando solo una matriz de valores de amplitud de la traza 1000

[Traza_1000, sampint, textheader] = altreadsegy('J:\Estudio\Maestria\Tesis\DATA\QUI-07-17_Traza_1000.segy','textheader','yes');

%Se realiza el llamado de nuevo de los segy para crear matriz de valores de amplitud de los segy de los gathers pre-apilados
[QUI_07_17, sampint17, textheader17] = altreadsegy('J:\Estudio\Maestria\Tesis\DATA\QUI-07-17_GATHERS.segy','textheader','yes');
[QUI_07_18, sampint18, textheader18] = altreadsegy('J:\Estudio\Maestria\Tesis\DATA\QUI-07-18_GATHERS.segy','textheader','yes');


%se crea una variable con tiempo desde 0 hasta 4 seg para poder graficar la tabla del dataout
t = 0:sampint:(size(Traza_1000,1)-1)*sampint;

%guardar la rata de muestreo
rm=sampint;

%se plotea el resultado
figure(1);
plot(t,Traza_1000);figure(gcf);
figure(2);
wiggle(Traza_1000);

%crea un archivo para guardar la matriz dela ondicula
savefile='Traza_1000';
save(savefile, 'Traza_1000','rm','t') % Guarda la ondicula y el rm

savefile='QUI_07_17';
save(savefile, 'QUI_07_17','rm','t') % Guarda la ondicula y el rm

savefile='QUI_07_18';
save(savefile, 'QUI_07_18') % Guarda la ondicula y el rm







