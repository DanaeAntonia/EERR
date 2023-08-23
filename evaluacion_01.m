%% cargo mis datos
clear all clc
a=readmatrix('horas.csv'); 
vientos=a(1:365*24,2)*(1000/(60*60)); %vientos en m/s %en horas
dir=a(1:365*24,3); %direccion del viento 
%% 2.a dibujar rosa de los vientos 
%%
%separo por meses.
dm=[0*24 31*24 28*24 31*24 30*24 31*24 30*24 31*24 31*24 30*24 31*24 30*24 31*24]; % el primero es cero por que tiene que partir el día 1
dms=cumsum(dm);
% vectores correspondientes a cada mes
ene = [vientos(1:dms(2)),dir(1:dms(2))]; 
feb = [vientos(dms(2)+1:dms(3)),dir(dms(2)+1:dms(3))]; 
mar = [vientos(dms(3)+1:dms(4)),dir(dms(3)+1:dms(4))]; 
abr = [vientos(dms(4)+1:dms(5)),dir(dms(4)+1:dms(5))];
may = [vientos(dms(5)+1:dms(6)),dir(dms(5)+1:dms(6))]; 
jun = [vientos(dms(6)+1:dms(7)),dir(dms(6)+1:dms(7))];
jul = [vientos(dms(7)+1:dms(8)),dir(dms(7)+1:dms(8))];
ago = [vientos(dms(8)+1:dms(9)),dir(dms(8)+1:dms(9))]; 
sep = [vientos(dms(9)+1:dms(10)),dir(dms(9)+1:dms(10))]; 
oct = [vientos(dms(10)+1:dms(11)),dir(dms(10)+1:dms(11))]; 
nov = [vientos(dms(11)+1:dms(12)),dir(dms(11)+1:dms(12))]; 
dic = [vientos(dms(12)+1:dms(13)),dir(dms(12)+1:dms(13))];
%realizo vector de octubre a marzo 
OaM_v=[oct(:,1) ; nov(:,1) ; dic(:,1);  ene(:,1); feb(:,1); mar(:,1)];
OaM_d=[oct(:,2) ; nov(:,2) ; dic(:,2);  ene(:,2); feb(:,2); mar(:,2)];
%realizo vector de abril a septiembre 
AaS_v=[abr(:,1) ; may(:,1); jun(:,1); jul(:,1); ago(:,1); sep(:,1)];
AaS_d=[abr(:,2) ; may(:,2); jun(:,2); jul(:,2); ago(:,2); sep(:,2)];
%% rosa de octubre a marzo 
figure()
wind_rose(OaM_d,OaM_v)
title('Rosa de los vientos estación Aeropuerto Carlos Ibañes, Punta Arenas, de Octubre a Marzo 2022')
set(gcf,'color','w')

%% rosa de abril a sep
figure()
wind_rose(AaS_d,AaS_v)
title('Rosa de los vientos estación Aeropuerto Carlos Ibañes, Punta Arenas, de Abril a Septiembre 2021')
set(gcf,'color','w')

%% Hacemos un histogrma con la magnitud de la velocidad del viento 
figure()
histogram(vientos,30, 'EdgeColor', 'black', 'FaceColor', 'blue')
xlabel('Magnitud de la velocidad del viento [m/s]')
ylabel('Frecuencia')
title('Histograma de la magnitud de la velocidad del viento, Punta Arenas')
set(gcf,'color','w')
%% Grafico tipo Chascón
%% 2.b Diagrama de vector progresivo
%% 
% descompones velocidad en vx y vy
%vx = -vo sin(thetha*pi/180)
%vy = -vo cos(thetha*pi/180)
%donde 
% vo = vh_2021
% theta = dh_2021
%paso a m/s ya que mis datos estan en km/h
vx= -vientos.*sin(dir.*pi/180)*3600; %multiplico por tiempo para obtener[m]
vy= -vientos.*cos(dir.*pi/180)*3600;
%sumo los elementos anteriores con función cumsum
x = cumsum(vx);
y = cumsum(vy);
%realizo diagrma
figure()
plot(x,y)
xlabel('Metros')
ylabel('Metros')
grid on 
xlim([-0.4*10^7 4*10^7])
ylim([-0.4*10^7 4*10^7])
set(gcf,'color','w')
title('Diagrama de vector progresivo de los vientos durante un año')

%% Calcule los parámetros necesarios y muestre gráficamente la distribución de Weibull para los datos sin normalizar y para los valores
%normalizados. Muestre sus cálculos. 
%% Nomalizado 
%ordenamos los datos
v_ord=sort(vientos);
indices_nan = isnan(v_ord);
% Reemplazar los elementos NaN por 0
v_ord(indices_nan) = [];

%normalizamos los datos ordenados 
x=v_ord./(mean(v_ord));
%calculamos k y c
%la k se calcula con los datos normalizados 
k = (std(x)./mean(x))^-1.086;
% c se calcula con la funcion gamma
c = 1/gamma(1+1/k);
%formula de la distribucion de weibull
p = (k/c).*(x/c).^(k-1).*exp(-(x./c).^k); 
figure()
plot(x,p,'LineWidth',3)
% hold on 
% histogram(x,30, 'EdgeColor', 'black', 'FaceColor', 'blue')
title('Distribución de Weibull para Punta Arenas 2022')
xlabel('Velocidad normalizada (v/|v|)')
ylabel('P(v/|v|)')
grid on
set(gcf,'color','w')
axis tight
%% Sin normalizar
% Ordenamos los datos
v_ord = sort(vientos);
indices_nan = isnan(v_ord);
% Reemplazar los elementos NaN por 0
v_ord(indices_nan) = [];
% Calculamos k y c
k = (std(v_ord) / mean(v_ord))^-1.086;
c = mean(v_ord) / gamma(1 + 1/k);
% Generamos un rango de valores para x
x = linspace(min(v_ord), max(v_ord), 100);
% Calculamos la distribución de Weibull para datos sin normalizar
p = (k / c) * (x / c).^(k-1) .* exp(- (x / c).^k);
% Realizamos la gráfica de la distribución
figure()
plot(x, p, 'LineWidth', 3)
title('Distribución de Weibull para Punta Arenas 2022')
xlabel('Velocidad de viento')
ylabel('P(v)')
grid on
set(gcf, 'color', 'w')
axis tight
%% Codigo entregado y el mío.
clear all
close all
clc
%% EXTRACT AND PLOT RAW DATA
a=readmatrix('horas.csv'); 
v=a(1:365*24,2)*(1000/(60*60));
indices_nan = isnan(v);
% Reemplazar los elementos NaN por 0
v(indices_nan) = [];
a=readmatrix('horas.csv'); 
vientos=a(1:365*24,2)*(1000/(60*60)); %vientos en m/s %en horas
%%
% Linearize distributions (see papers)
ln = log(uniqueVals);
lnln = log(-log(1-cumFreq));

% Check wether the vectors contain inifinite values, if so, remove them
test = isinf(lnln);
for i=1:nbUniqueVals
    if (test(i)==1)
        ln(i)= [];
        lnln(i)= [];
    end
end

% Extract the line parameters (y=ax+b) using the polyfit function
params = polyfit(ln,lnln',1);
a = params(1);
b = params(2);
y=a*ln+b;

% Compare the linealized curve and its fitted line
figure
plot(ln,y,'b',ln,lnln,'r')
title('Linearized curve and fitted line comparison');
xlabel('x = ln(v)');
ylabel('y = ln(-ln(1-cumFreq(v)))');
set(gcf,'color','w')
% Extract the Weibull parameters c and k
k = a
c = exp(-b/a)
%% CHECK RESULTS
% Define the cumulative Weibull probability density function
% F(V) = 1-exp(-((v/c)^k)) = 1-exp(-a2), with a1 = v/c, a2 = (v/c)^k
a1 = uniqueVals/c;
a2 = a1.^k;
cumDensityFunc = 1-exp(-a2); 
% Define the Weibull probability density function
%f(v)=k/c*(v/c)^(k-1)*exp(-((v/c)^k))=k2*a3.*exp(-a2), 
% with  k2 = k/c, a3 = (v/c)^(k-1)
k1 = k-1;
a3 = a1.^k1;
k2 = k/c;
densityFunc = k2*a3.*exp(-a2);  
% Plot and compare the obtained Weibull distribution with the frequency plot
% Ordenamos los datos
v_ord = sort(vientos);
indices_nan = isnan(v_ord);
% Reemplazar los elementos NaN por 0
v_ord(indices_nan) = [];
% Calculamos k y c
k_sn = (std(v_ord) / mean(v_ord))^-1.086;
c_sn = mean(v_ord) / gamma(1 + 1/k_sn);
% Generamos un rango de valores para x
x = linspace(min(v_ord), max(v_ord), 100);
% Calculamos la distribución de Weibull para datos sin normalizar
p = (k_sn / c_sn) * (x / c).^(k_sn-1) .* exp(- (x / c_sn).^k_sn);
figure()
subplot 211
plot(uniqueVals,prob,'k.',uniqueVals,densityFunc, 'r','LineWidth',1)
hold on 
plot(x, p,'b', 'LineWidth', 1)
title('Weibull probability density function');
xlabel('v');
ylabel('f(v)');
legend('Datos','Weibull enfoque gráfico','Weibull enfoque std,mean y gamma')
set(gcf,'color','w')
%% Actualizar los Gráficos. 
%% leo datos
a=xlsread('CEN-hist_gen_de_energia_por_tecnologia.xlsx',2); %coordinadora ...
%generación de energía por tecnología
anos=a(:,1)';
e=a(:,2:14)';
hidr=a(:,2); %primera
carbon=a(:,3); %segunda
Diesel=a(:,4); %tercera
GasNat=a(:,5); %cuarta
oil=a(:,6);%octavae
petcoke=a(:,7);%sexta
cogeneracion=a(:,8);
biogas=a(:,9); %novena
biomasa=a(:,10);%quinta
Eolica=a(:,11);%septime
solar=a(:,12);%decima
termosolar=a(:,13);%doceava
geotermica=a(:,14);%onceava
%% ordenando energías 
energias=[hidr,cogeneracion,biomasa,Eolica,biogas,solar,geotermica,termosolar,carbon,Diesel,GasNat,petcoke,oil];
%% caluclo porcentajes 
for i= 1:23 %año 2000-2022
  total=sum(energias(i,:));
  ren=sum(energias(i,1:8));
  renovables(i)=ren
  porcentaje(i)=round((ren/total)*100,2);
end
%% grafico
colors = {'b', 'g',[0.4660, 0.6740, 0.1880] ,[0.3010, 0.7450, 0.9330] ,[0 0.5 0] , 'y',[1 0.5 0] ,'r', 'k',[0.2 0.2 0.2] ,[0.7, 0.7, 0.7] ,[0.5, 0.5, 0.5] ,[30/255 30/255 30/255] };
figure()
b= bar(anos', [energias(:,1:8), energias(:,9:end)], 'stacked');
for i = 1:numel(b)
    b(i).FaceColor = colors{i};
end
hold on
plot(anos, renovables, 'm-', 'LineWidth', 3)
hold off
legend('Hidráulica', 'Cogeneración', 'Biomasa', 'Eólica', 'Biogás', 'Solar', 'Geotérmica', 'Termosolar', 'Carbón', 'Diesel', 'Gas Natural','Petcoke', 'Petróleo','Total E.Renovables');
set(gca,'XTick',anos) % establece los valores de años como etiquetas del eje x
set(gca,'XTickLabelRotation',90) %
xlim([anos(1), anos(end)])
xlabel('Años','Fontsize',15)
title('Producción de energía por tipo de tecnología entro los años 2000 y 2022','Fontsize',20)
ylabel('Generación de energía (GWh)','Fontsize',15)
grid on
set(gcf,'color','w')
%% lo mismo de lo anterior pero en porcentajes 
for i=1:23 
    for j=1:13
        total=sum(energias(i,:));
        porcentaje(i,j)= round((energias(i,j)/total)*100,2);
        ren=sum(energias(i,1:8));
        renovables(i)=round((ren/total)*100,2);
        
    end 
end
%% grafico en porcentajes
colors = {'b', 'g', [0.4660, 0.6740, 0.1880], [0.3010, 0.7450, 0.9330], [0 0.5 0], 'y', [1 0.5 0], 'r', 'k', [0.2 0.2 0.2], [0.7, 0.7, 0.7], [0.5, 0.5, 0.5], [30/255 30/255 30/255]};
figure()
b = bar(anos', [porcentaje(:,1:8), porcentaje(:,9:end)], 'stacked');
for i = 1:numel(b)
    b(i).FaceColor = colors{i};
end
hold on
plot(anos, renovables, 'm-', 'LineWidth', 3)
legend('Hidráulica', 'Cogeneración', 'Biomasa', 'Eólica', 'Biogás', 'Solar', 'Geotérmica', 'Termosolar', 'Carbón', 'Diesel', 'Gas Natural', 'Petcoke', 'Petróleo', 'Total E.Renovables', 'FontSize', 5);
% Agregar etiquetas de porcentaje
for i = 1:numel(anos)
    text(anos(i), renovables(i)+7, [num2str(renovables(i)) '%'], 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'Rotation', 90, 'FontSize', 12, 'Color', 'white');
end
set(gca, 'XTick', anos) % establece los valores de años como etiquetas del eje x
set(gca, 'XTickLabelRotation', 90)
xlim([anos(1), anos(end)])
xlabel('Años', 'FontSize', 15)
title('Porcentaje de energía por tipo de tecnología entre los años 2000 y 2022', 'FontSize', 20)
ylabel('Porcentaje de energía', 'FontSize', 15)
grid on
set(gcf, 'color', 'w')
ylim([0, 100]);
ytickformat('%.2f%%');
set(gca, 'yaxis', struct('TickLabelFormat', '%.0f'))
