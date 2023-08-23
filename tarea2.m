%% cargo mis datos
a=readmatrix('Dichato_2021_horas.csv'); 
vientos=a(1:365*24,2); %vientos en m/s %en horas
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
%% 2.b Diagrama de vector progresivo
%% 
% descompones velocidad en vx y vy
%vx = -vo sin(thetha*pi/180)
%vy = -vo cos(thetha*pi/180)
%donde 
% vo = vh_2021
% theta = dh_2021
vientos.*(1000/3600); %paso a m/s ya que mis datos estan en km/h
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
xlim([-0.5*10^7 5*10^7])
ylim([-0.5*10^7 5*10^7])
set(gcf,'color','w')
title('Diagrama de vector progresivo de los vientos durante un año')

%% 2.c Muestre graficamente la distribucion de Weibull para esos datos.
%ordenamos los datos
v_ord=sort(vientos);
%normalizamos los datos ordenados 
x=v_ord./(mean(v_ord));
%calculamos k y c
%la k se calcula con los datos normalizados 
k = (std(x)./mean(x))^-1.086;
% c se calcula con la funcion gamma
c = 1/gamma(1+1/k);
%formula de la distribucion de weibull
p = (k/c).*(x/c).^(k-1).*exp(-(x./c).^k); 
%realizamos la grafica de la distriución 
figure()
plot(x,p,'LineWidth',3)
title('Distribucipon de Weibull para Dichato 2021')
xlabel('Velocidad normalizada (v/|v|)')
ylabel('P(v/|v|)')
grid on
set(gcf,'color','w')
axis tight
%% 2.d ¿Qué porcentaje del tiempo la intensidad del viento es mayor a 3 m/s?
t_d= 1-exp(-(3/mean(vientos))*gamma(1+(1/k)))^k;
%este t esta en cifra decimal, hay que multiplicarlo por 100
tp_d=t_d*100;
%% 2.e ¿Qué porcentaje del tiempo la intensidad del viento está entre 3 y 12 m/s?
t_e = (exp(-(3/mean(vientos))*gamma(1+(1/k)))^k) - (exp(-(12/mean(vientos))*gamma(1+(1/k)))^k);
%este t2 esta en decimal, hay que multiplicarlo por 100
tp_e=t_e*100;
%% 3.a Calcule los promedios mensuales de la intensidad de los vientos
%usados en el ejercicio anterior
%mis datos estan en horas por lo que debo calcular el promedio diario
n=0;
for i=1:365 %dias
    vd(i)=mean(vientos(n+1:n+24)); %promedio 
    n=n+24;
end
% ahora calulo los promedios mensuales
dm=[0 31 28 31 30 31 30 31 31 30 31 30 31]; %dias del mes 
%y empiezo del cero pq tiene que empezar del dato 1
dms=cumsum(dm); %se va sumando los dias del mes 
%esta suma indica la posicion del mes que viene
for i=1:12
    vm(i)=mean(vd(dms(i)+1:dms(i+1))); %promedios mensuales
end
meses = {'Ene', 'Feb', 'Mar', 'Abr', 'May', 'Jun', 'Jul', 'Ago', 'Sep', 'Oct', 'Nov', 'Dic'};
figure()
plot(vm)
hold on
plot(vm,'r.','MarkerSize',25)
xticks(1:length(meses));
xticklabels(meses);
title('Intensiad promedio del viento por mes');
ylabel('Intensidad promedio (m/s)');
grid on
set(gcf,'color','w')
axis tight
%%  haga una esquema de cómo varían los vientos con la altura
%esos vientos con la altura, hasta 80 m de altura, si la medici´on hubiese
%sido hecha en un potrero con pasto corto, y si hubiese sido hecha en un
%bosque.
%segun gpt esta a 10 metros de altura. 
z1=10; % metros altura a la que se tomaron los datos 
%z2=[35,36,37,38,39,40,41,42,43,44,45,46,47,48,49, 50,51,52,53,54, 55,56,57,58,59 60,61,62,63,64, 65,66,67,68,69, 70,71,72,73,74,75,76,77,78,79,80]; % altura a la cual quiero saber la velocidad
z2=[10:2:80];
%comenzamos en 35 ya que el largo de la aspa es de 35 metros
z0_p= 0.03; %pasto
z0_b= 1 ; %bosque
%vm=sum(vm)/length(vm); %promedio anual de la velocidad.
%formula % v(z2)/v(z1) = log(z2/z0) / log(z1/z0);
%v(z1) es la velocidad media por mes(vm)
for j=1:12
for i=1:length(z2)
v_p(i,j)=(log(z2(i)/z0_p)/log(z1/z0_p))*vm(j); % velocidad suelo con pasto
v_b(i,j)=(log(z2(i)/z0_b) / log(z1/z0_b))*vm(j); %velocidad suelo con bosque
end 
end
%en mi matriz de velocidad las columnas representan los meses  las filas as
%distintas alturas hasta 80 metros  
%% grafica para bosque
figure();
hold on;
for i=1:12
    if i==1 % Enero en azul
        plot(v_b(:,i),z2,'.b','MarkerSize',8);
    elseif i==2 % Febrero en rojo
        plot(v_b(:,i),z2,'.r','MarkerSize',8);
    elseif i==3 % Marzo en verde
        plot(v_b(:,i),z2,'.g','MarkerSize',8);
    elseif i==4 % Abril en magenta
        plot(v_b(:,i),z2,'.m','MarkerSize',8);
    elseif i==5 % Mayo en cian
        plot(v_b(:,i),z2,'.c','MarkerSize',8);
    elseif i==6 % Junio en amarillo
        plot(v_b(:,i),z2,'.y','MarkerSize',8);
    elseif i==7 % Julio en negro
        plot(v_b(:,i),z2,'.k','MarkerSize',8);
    elseif i==8 % Agosto en cyan oscuro
        plot(v_b(:,i),z2,'.','Color',[0 0.5 0.5],'MarkerSize',8);
    elseif i==9 % Septiembre en azul oscuro
        plot(v_b(:,i),z2,'.','Color',[0 0 0.5],'MarkerSize',8);
    elseif i==10 % Octubre en verde oscuro
        plot(v_b(:,i),z2,'.','Color',[0 0.5 0],'MarkerSize',8);
    elseif i==11 % Noviembre en rojo oscuro
        plot(v_b(:,i),z2,'.','Color',[0.5 0 0],'MarkerSize',8);
    else % Diciembre en morado oscuro
        plot(v_b(:,i),z2,'.','Color',[0.5 0 0.5],'MarkerSize',8);
    end
end
xlabel('Velocidad del viento (m/s)');
ylabel('Altura (m)');
title('Velocidad del viento mensual a diferentes alturas, bosque');
legend(meses);
grid on 
set(gcf,'color','w')
axis tight
%% Ahora para pasto corto
figure();
hold on;
for i=1:12
    if i==1 % Enero en azul
        plot(v_p(:,i),z2,'.b','MarkerSize',8);
    elseif i==2 % Febrero en rojo
        plot(v_p(:,i),z2,'.r','MarkerSize',8);
    elseif i==3 % Marzo en verde
        plot(v_p(:,i),z2,'.g','MarkerSize',8);
    elseif i==4 % Abril en magenta
        plot(v_p(:,i),z2,'.m','MarkerSize',8);
    elseif i==5 % Mayo en cian
        plot(v_p(:,i),z2,'.c','MarkerSize',8);
    elseif i==6 % Junio en amarillo
        plot(v_p(:,i),z2,'.y','MarkerSize',8);
    elseif i==7 % Julio en negro
        plot(v_p(:,i),z2,'.k','MarkerSize',8);
    elseif i==8 % Agosto en cyan oscuro
        plot(v_p(:,i),z2,'.','Color',[0 0.5 0.5],'MarkerSize',8);
    elseif i==9 % Septiembre en azul oscuro
        plot(v_p(:,i),z2,'.','Color',[0 0 0.5],'MarkerSize',8);
    elseif i==10 % Octubre en verde oscuro
        plot(v_p(:,i),z2,'.','Color',[0 0.5 0],'MarkerSize',8);
    elseif i==11 % Noviembre en rojo oscuro
        plot(v_p(:,i),z2,'.','Color',[0.5 0 0],'MarkerSize',8);
    else % Diciembre en morado oscuro
        plot(v_p(:,i),z2,'.','Color',[0.5 0 0.5],'MarkerSize',8);
    end
end
xlabel('Velocidad del viento (m/s)');
ylabel('Altura (m)');
title('Velocidad del viento mensual a diferentes alturas, pasto corto');
legend(meses);
grid on 
set(gcf,'color','w')
axis tight
%% si tomo el promedio anual de vm y calculo para diferentes alturas en los dos casos
vm=sum(vm)/length(vm); %promedio anual 

for i=1:length(z2)
va_p(i)=(log(z2(i)/z0_p)/log(z1/z0_p))*vm; % velocidad suelo con pasto
va_b(i)=(log(z2(i)/z0_b) / log(z1/z0_b))*vm; %velocidad suelo con bosque
end 
figure()
plot(va_b,z2,'.r','MarkerSize',9)
hold on
plot(va_p,z2,'.b','MarkerSize',9)
xlabel('Intesidad del viento (m/s)')
ylabel('Altura (m)')
legend('Bosque','Pasto corto')
grid on
title('Velocidad del viento anual a diferentes alturas')
set(gcf,'color','w')
axis tight
%% 3.b caluclar potencia y pasar a dinero. 
%Voy a usar valor promedio anual por altura usare va_p y va_b
%entonces vamos  suponer que las aspas comienzan a la altura 10 hasta 80 
%es decir el diametro del aerogenerador es de 70 metros desde 10 a 80 es 
%decir a una altura de 45 metros esta instalado el aerogenerador.
%por eso tengo alturas desde 10:80 cada 2 metros por lo tanto vamos a
%suponer que el area que barre el aerogenerador es un cadrado de 70x70
%pero ese cuadrado lo vamos a dividr en 36 ya que son todas las alturas que
%tengo entonces para calcular los potenciales tendremos 36 potenciales en
%36 areas, es decir el cuadrado lo parto en 36 partes de 70x2 para cada
%altura (36) de 10 metros hasta 80 metros.
radio=35; %m
area=2*radio*2; %m^2
rho=1.29;
cp=0.4; %el gamma me da la eficiencia pero igual pregunta al profe 
%despues sumo todas mis potencias y calculo mi valor de potencia 
%según boleta de luz KWh esta en 127.4 pesos entonces el MWh esta en
%127.400 pesos
valor= 127400.% valor del MWh 
for i=1:length(va_b)
    M_b(i)=(((1/2)*rho*area*va_b(i)^3*cp)) / 1000000; %MW
    M_p(i)=(((1/2)*rho*area*va_p(i)^3*cp)) / 1000000
end
potencia_b=sum(M_b)*24*365 %MWh %encontramos potencia promedio
money_b=(potencia_b)*valor %MWh
potencia_p=sum(M_p)*24*365 %MWh %encontramos potencia promedio
money_p=(potencia_p)*valor
%% otra forma 
%vamos a calcular una velocidad promedio y con eso calcularemos la potencia
vm_b=sum(va_b)/length(va_b);
vm_p=sum(va_p)/length(va_p);
Mm_b=((((1/2)*rho*area*vm_b^3*cp))*24*365) / 1000000; %MWh
Mm_p=((((1/2)*rho*area*vm_p^3*cp))*24*365) / 1000000; %MWh
money_mb=Mm_b*valor
money_mp=Mm_p*valor

 
 
 
 
 
%debo calcular promedio anual pafra cada altura
for j=1:length(z2)
        vm_p(j)=sum(v_p(j,1:12))/12; %promedia anual a diferentes alturas
        vm_b(j)=sum(v_b(j,1:12))/12;
end
figure()
plot(vm_b,z2,'.r','MarkerSize',9)
hold on
plot(vm_p,z2,'.b','MarkerSize',9)
xlabel('Intesidad del viento (m/s)')
ylabel('Altura (m)')
legend('Bosque','Pasto corto')
grid on
title('Velocidad del viento anual a diferentes alturas')
set(gcf,'color','w')
axis tight
%ahora que tendo las velocidades anuales para cada altura calculo la
%potencia de un aerogenerador de tres aspas de 35 metros
radio=35; %m
area=pi*radio^2; %m^2
rho=1.29;
cp=0.43;
valor=127.4; %clp
for i=1:length(vm_b)
    M(i)=(((1/2)*rho*area*cp*(gamma(1+(3/k))/(gamma(1 + (1/k)))^3)*vm_b(i)^3*0.43)/1000)*24*365; %KWh
    money(i)=(M(i)*valor)/1000; %MWh
end
% el valor del kwh en chile es de 127.4 según boleta de luz de mi casa 
%% FIN
%% ME FALTA EL EJERCICIO DISEÑADO 
%Compare la rosa de los vientos de una estación al norte y otra en el sur
%para el mes de enero y mayo.








%% 3.b calcular la potencia en cada caso 
% dato radio= 35 metros
%potencia no perturabda es M=(1/2)*rho*area*viento^3 
%area=pi*radio^2 
radio=35; %m
area=pi*radio^2; %m^2
rho=1.29; %desidad del aire en kg/m3
cd=0.07;
M=(((1/2)*rho*area*cd*(gamma(1+(3/k))/(gamma(1 + (1/k)))^3)*vm^3)/1000000)*24*365; %vm esta en m/s 

% esto esta en MW promedio en un mes para saber la energia electrica
% generada por hora es decir MWh debo saber las horas del mes y multiplicar
% por esas horas 
dm=[31*24 28*24 31*24 30*24 31*24 30*24 31*24 31*24 30*24 31*24 30*24 31*24]; % el primero es cero por que tiene que partir el día 1
dms=cumsum(dm);
for i=1:12
    MWh(i)= M(i)*dm(i);
end
% pasamos a dinero considerando que 1MWh vale Según datos de la Comisión...
%Nacional de Energía de Chile, el precio promedio ponderado de la energía...
%eléctrica en el mercado spot (mercado de corto plazo) para el año 2021...
%fue de aproximadamente 44,2 dólares estadounidenses por MWh.
% suponemos que el dolar esta a 812 pesos 44,2*812 1 MWh en pesos chilenos
dinero=sum(MWh*44.2*812);
%%
%paso de hora a día
%n=0;
%for i=1:365
%    vientos_d=mean(vientos(n+1:n+24,1)); %valores diarios
%    n=n+24;
%end
%del vector de dias extraer los meses de maarzo mes 3 (0 31 28)
%for i= 3:10 
%    v_MaO=vientos_d(dms(i):dms(i+1));
%end
%pasar de dias a meses
dm=[0 31 28 31 30 31 30 31 31 30 31 30 31]; % el primero es cero por que tiene que partir el día 1
dms=cumsum(dm);
vientos_MaO= vientos(60:273,1); %en la posición 60 se encuentra 1 de marzo y en 273 el ultimo dia de octubre
dir_MaO=dir(60:273,1);
%% realizamos rosa de los vientos para periodo de marzo a octubre
wind_rose(dir_MaO,vientos_MaO)

for i=1:12
    vientos_m(i)=vientos(dms(i):dms(i+1),1); %valores mensuales
end

%% conocer velocidad en superficie
 %si conoces el viento a cierta altura podemos estimar la velocidad del
 %viento paralelo a la superficie 
%z1 altura a la que se toman los datos
%z2 a la altura a la que quiero conocer la velocidad del viento
%z0 si es un prado o un bosque, el valor debemos mirarlo en las clases del
%profe 
% V(z2)/V(z1)= log(z2/z0)/log(z1/z0) 
%% potencia 
% uno tiene que dar la eficiencia, buscar en internet
%calcular energia cinetica de areogenerador 1/2(masa)'(velocidad)^2
%masa' m'=rho*volumen', volumen'=rho*area*velocidad 
%metemos masa en ecuación de energía cinetica y obtenemos E_c=1/2
%*rho*area*velocidad^3 
%pero eso no es todo debemos agregar Función Gamma y tendremnos
%Ec=1/2 rho*area*velocidad^3* gamma(1+3/k)/gamma(1+1/k)^3 * eficiencia del
%aerogenerador   
% potencia en watt a kw= potencia/1000, MW= potencia /1000000;
% kwhora*cuanto vale kilowat hora= Dinero.