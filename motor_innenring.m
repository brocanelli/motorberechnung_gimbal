% Berechnung für Antriebsmoment
% Trapezgewinde
close all;
clear;
clc;

kv = 355; %rpm/V % Masse 275g
P = 2*1e-3; % Steigung in mm

%% Kinematik
d = 64.5*2*1e-3; %Angriffspunkt des Gimbal-Motors
x0 = 65*1e-3; %Abstand des oberen Gabelpunktes zum Globaeln Koordinatensys.
l = 40*1e-3; %Länge der Gabel
a = 19*1e-3; %Höhe unteres gabelauge

syms phi
% phi = linspace(-10,10,21)*pi/180;
psi = asin(x0/l - cos(phi)*d/2/l - a/l*sin(phi));
f_psi = matlabFunction(psi);

% gamma = [-10 -9 -8 -7 -6 -5 -4 -3 -2 -1 0 1 2 3 4 5 6 7 8 9 10];
% alpha_innen = [6.862 6.122 5.411 4.727 4.070 3.442 2.841 2.268 1.722 1.205 0.716... 90ig
%     0.255 -0.177 -0.581 -0.957 -1.304 -1.623 -1.912 -2.173 -2.404 -2.607]*(-1);
% % alpha_ausen = [

%% Dynamik
Fs = 10e3; %Schubkraft
c = 20*1e-3; %Exzentrizität
Isd = 254675.319*1e-6; %kg mm^2
Isi = 2363.620*1e-6; %kg mm^2 (Relativer Fehler = 0.001992%)
Is = Isd + Isi;
A = 10*pi/180;
f = 2; %Hz
t = linspace(0,1,101);

phi_t = A*sin(2*pi*f*t);
phi_t_dot = 2*pi*f*A*cos(2*pi*f*t);
phi_t_ddot = -4*pi^2*f^2*A*sin(2*pi*f*t);

Fst1 = -2*(Is*phi_t_ddot + Fs*c)./(cos(f_psi(phi_t)).*cos(phi_t)+sin(f_psi(phi_t)).*sin(phi_t))/d;
Fst2 = -2*(Is*phi_t_ddot - Fs*c)./(cos(f_psi(phi_t)).*cos(phi_t)+sin(f_psi(phi_t)).*sin(phi_t))/d;

%% Berechnung
d_1 = 10*1e-3; %Durchmesser des Trapezgewindes in mm
H1 = 0.5*P; % Tiefe des Trapezgewindes in mm
beta = 30*pi/180; % Flankenwinkel

%% Ausrechnen der Geschwindigkeit auf y-Achse
% Umrechnen auf Winkelgeschwindigkeit
psi_dot = (d/2*sin(phi_t)-a*cos(phi_t)).*phi_t_dot./l./cos(f_psi(phi_t));
vky = -d/2*cos(phi_t).*phi_t_dot - a*sin(phi_t).*phi_t_dot - l*sin(phi_t).*psi_dot;
gamma_dot = 2*pi*vky/P;
n = gamma_dot/2/pi; %rps

%% --------
% mu_g = 0.35; %worst case Reibung Gewinde
% mu_l = 0.35; %worst case Reibung Lager
mu_g = 0.09; %worst case Reibung Gewinde
mu_l = 0.09; %worst case Reibung Lager

mu_g = 0.01; %worst case Reibung Gewinde
mu_l = 0.01; %worst case Reibung Lager

d2 = d_1-H1; %Flankendruchemesser
r2 = d2/2; % Flankenradius
% rL = d2/2; % Lagerdurchmesser (dw+dh)/4
rL = (11+7)/4*(1e-3); %Roloff-Matek

phi = atan(P/d2/pi);
rho = atan(mu_g/cos(beta/2));

% ------------------------------------------------------------------------
Fay1 = Fst1.*cos(f_psi(phi_t));
Fay2 = Fst2.*cos(f_psi(phi_t));
Fax1 = Fst1.*sin(f_psi(phi_t));
Fax2 = Fst2.*sin(f_psi(phi_t));

% MG1 = Fay1*tan(phi-rho)*r2; % Gewindereibmoment
% ML1 = Fay1*mu_l*rL;
% MA1 = abs(MG1) + abs(ML1);

% MG1 = Fay1*abs(tan(phi-rho))*r2; % Gewindereibmoment
% ML1 = Fay1*mu_l*rL;
% MA1 = MG1 + ML1;
%
% MG2 = Fay2*abs(tan(phi-rho))*r2; % Gewindereibmoment
% ML2 = Fay2*mu_l*rL;
% MA2 = MG2 + ML2;
if phi-rho < 0
    sprintf('ACHTUNG: SELBSTHEMMUNG')
end

MA1 = [];
MA2 = [];

for idx = 1:length(phi_t) % falle 1
    if phi_t_dot(idx) > 0 % Senken
%         fprintf('senken \n');
        MA1 = [MA1; Fay1(idx)*(tan(phi-rho))*r2 + Fay1(idx)*mu_l*rL];
        MA2 = [MA2; Fay2(idx)*(tan(phi-rho))*r2 + Fay2(idx)*mu_l*rL];
    else % heben
        MA1 = [MA1; Fay1(idx)*(tan(phi+rho))*r2 + Fay1(idx)*mu_l*rL];
        MA2 = [MA2; Fay2(idx)*(tan(phi+rho))*r2 + Fay2(idx)*mu_l*rL];
    end
end

% Gewindereibmoment nach Rolof-Matek
% Lagerreibmoment nach FAG
% AXK0821-TV
% https://medias.schaeffler.com/medias/de!hp.tg.cat/tg_hr*ST4_18687224075
% f0 = 3;
% f1 = 0.0015;
% dm = (21+8)/2; %Mittlerer Lagerdurchmesser
% mu = 100; %mm^2/s kinematische Viskosität bei 40°C ISOVG100
% 
% for idx = 1:length(phi_t) % falle 1
%     if abs(mu*n*60) >= 2000
%         M0 = f0*(mu*n(idx)*60)^(2/3)*dm^3*10^(-7)*1e-3;
%     else
%         M0 = f0*160*dm^3*10^(-7)*1e-3;
%     end
%     
%     M1 = f1*Fay1(idx)*dm*1e-3; %Nm
%     MR = M1+M0;
%     
%     if phi_t_dot(idx) > 0 % Senken
% %         fprintf('senken \n');
%         MA1 = [MA1; Fay1(idx)*(tan(phi-rho))*r2 + MR];
%         MA2 = [MA2; Fay2(idx)*(tan(phi-rho))*r2 + Fay2(idx)*mu_l*rL];
%     else % heben
%         MA1 = [MA1; Fay1(idx)*(tan(phi+rho))*r2 + MR];
%         MA2 = [MA2; Fay2(idx)*(tan(phi+rho))*r2 + Fay2(idx)*mu_l*rL];
%     end
% end

%% Calculation of required Power and current draw
Pow1 = gamma_dot.*MA1'; % Watt
Pow2 = gamma_dot.*MA2'; % Watt

%% Maximum Current
V = n*60/kv;
I1 = MA1*2*pi/60*kv;
I2 = MA2*2*pi/60*kv;

%% Plot figures ----------------------------------------------------------
close all;

figure(1);
subplot(5,1,1)
% title('moment over angle');
yyaxis left;
plot(t,MA1);hold on;
yyaxis right;
plot(t,MA2);
% ylim([1.4,1.6]);
legend('MA1 [Nm]','MA2 [Nm]','Location','northwest');
hold off;

subplot(5,1,2)
% title('force over angle');
yyaxis left;
% ylim([3000,3400]);
plot(t,Fay1);hold on;
yyaxis right;
% ylim([-3400,-3000]);
plot(t,Fay2);
legend('FA1 [N]','FA2 [N]','Location','northwest');
hold off;

subplot(5,1,3)
% title('force over angle');
yyaxis left;
% ylim([3000,3400]);
plot(t,n*60);hold on;
yyaxis right;
% ylim([-3400,-3000]);
plot(t,V);
legend('n [rpm]','U [V]','Location','northwest');
hold off;

subplot(5,1,4)
% title('force over angle');
yyaxis left;
% ylim([3000,3400]);
plot(t,Pow1);hold on;
yyaxis right;
% ylim([-3400,-3000]);
plot(t,Pow2);
legend('P1 [W]','P2 [W]','Location','northwest');
hold off;

subplot(5,1,5)
% title('force over angle');
yyaxis left;
% ylim([3000,3400]);
plot(t,I1);hold on;
yyaxis right;
% ylim([-3400,-3000]);
plot(t,I2);
legend('I1 [A]','I2 [A]','Location','northwest');
hold off;

%% Output ----------------------------------------------------------------
stringM = ["------------------------------------";
    "\n";
    "Maximales Drehmoment: %0.5f Ncm";
    "\n";
    "------------------------------------"];

stringF = ["------------------------------------";
    "\n";
    "Maximale Kraft: %0.5f N";
    "\n";
    "------------------------------------"];

stringP = ["------------------------------------";
    "\n";
    "Maximale Leistung: %0.5f W";
    "\n";
    "------------------------------------"];

stringI = ["------------------------------------";
    "\n";
    "Maximaler Strom: %0.5f A";
    "\n";
    "------------------------------------"];

sprintf(join(stringM), max([max(max(MA1)),max(abs(MA2))])*1e2 )
sprintf(join(stringF), max([max(abs(Fay1)),max(abs(Fay2))]))
sprintf(join(stringP), max([max(abs(Pow1)),max(abs(Pow2))]))
sprintf(join(stringI), max([max(abs(I1)),max(abs(I2))]))


