load('EC_series (1).mat');
tiempo= datetime(time, 'ConvertFrom', 'datenum');
%% tomamos solo datos del año 2017 
tiempo_2012=tiempo(12648:13002);
value_2012=value(12648:13002);
%% a) caudal diria y clasificado
%% graficamos caudales diarios para año 2017
figure()
plot(tiempo_2012,value_2012,'b','LineWidth',3)
xlabel('Tiempo')
ylabel('Caudal m3/s')
title('Caudal diario del Río Bío-Bío desembocadura, año 2012')
grid on 
set(gcf,'color','w')
axis tight
%% Graficamos caudal clasificado 
clasificado = sort(value_2012, 'descend');
figure()
plot(clasificado,'b','LineWidth',3)
xlabel('Número de días')
ylabel('Caudal m3/s')
title('Caudal clasificado del Río Bío-Bío desembocadura, año 2012')
grid on 
set(gcf,'color','w')
axis tight
%% b) 
%% obtenemos caudal de equipamiento 
%calculamos el caudal ecologico, se lo restamos al caudal clasificado y
%luego obtenemos un caudal de equipamiento entre 80-100 días 
ecologico=mean(value)*0.1;
resta=clasificado-ecologico;
equipamiento=resta(90);
figure()
plot(clasificado,'b','LineWidth',3)
hold on 
plot(resta,'r','LineWidth',3)
plot(90,equipamiento,'ok','Linewidth',7)
xlabel('Número de días')
ylabel('Caudal m3/s')
title('Caudal de equipamiento del Río Bío-Bío desembocadura, año 2012')
grid on 
set(gcf,'color','w')
axis tight
legend('Clasificado','Clasificado-Ecologico','Equipamiento')
%% Ejercicio propuesto
%% leo datos
a=xlsread('CEN-hist_gen_de_energia_por_tecnologia.xlsx',2); %coordinadora ...
%generación de energía por tecnología
anos=a(:,1)';
e=a(:,2:14)';
hidr=a(:,2); %primera
figure()
bar(anos,hidr)
xlabel('Años')
ylabel('Generación de energía (GWh')
title('Energía hidromotriz en Chile')
grid on 
set(gcf,'color','w')
axis tight


