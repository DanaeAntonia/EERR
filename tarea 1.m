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

%% Ahora en barras separadas
colors = {'b', 'g',[0.4660, 0.6740, 0.1880] ,[0.3010, 0.7450, 0.9330] ,...
    [0 0.5 0] , 'y',[1 0.5 0] ,'r', 'k',[0.2 0.2 0.2] ,[0.7, 0.7, 0.7] ,...
    [0.5, 0.5, 0.5] ,[30/255 30/255 30/255] };
figure()
b= bar(anos', energias,5);
for i = 1:numel(b)
    b(i).FaceColor = colors{i};
end
hold off
legend('Hidráulica', 'Cogeneración', 'Biomasa', 'Eólica', 'Biogás', ...
    'Solar', 'Geotérmica', 'Termosolar', 'Carbón', 'Diesel', 'Gas Natural',...
    'Petcoke', 'Petróleo','Fontsize',6);
set(gca,'XTick',anos) % establece los valores de años como etiquetas del eje x
set(gca,'XTickLabelRotation',90) %
xlim([anos(1), anos(end)])
xlabel('Años','Fontsize',15)
title('Producción de energía por tipo de tecnología entre los años 2000 y 2022','Fontsize',20)
ylabel('Generación de energía (GWh)','Fontsize',15)
grid on
set(gcf,'color','w')
%% https://datosmacro.expansion.com/energia-y-medio-ambiente/emisiones-co2/chile
a=readtable('emicionesCo2.txt');
datos=table2array(emicionesCo2);

figure()
subplot 211
scatter(datos(:,1),datos(:,4))
ylabel('C02 (t)')
xlabel('Años')
title('Emiciones de C02 per capita en Chile')
grid on
set(gcf,'color','w')
xlim([datos(end,1) , datos(1,1)])
subplot 212
scatter(datos(:,1),datos(:,2), 'r')
ylabel ('C02 (Mt)')
xlabel('Años')
title('Emiciones de C02 por año en Chile')
xlim([datos(end,1) , datos(1,1)])
grid on
%% https://datos.bancomundial.org/indicator/EN.ATM.CO2E.KT
a=xlsread('API_EN.ATM.CO2E.KT_DS2_es_excel_v2_5348729.xls',1);
load('APIEN.mat')
anos=table2array(APIEN);
a2019=a(:,30); %emisiones de C02 para 265 paises en el año 2019
nombres(1,1)=[];


figure()
x = categorical(nombres)
bar(x,a2019)
bar(a2019)
set(gca,'xticklabel',nombres)
