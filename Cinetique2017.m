function  [DPn,C,ri] = Cinetique2017(a,b,c,d,e)


%%
%Cette section contient toutes les options de simulations choisies, 0 pour non, 1
%pour oui. On a, dans l'ordre, l'effet Trommsdorf, l'effet de
%vitrification, l'ajout d'un agent de transfert, l'épuisement de l'amorçeur
%et l'impact de la température. Ce variables sont utilisées plus tard dans
%le programme

global trommsdorf vitrification transfert ricst temperature
trommsdorf    = a;
vitrification = b;
transfert     = c;
ricst         = d;
temperature   = e;

%%
%Choix d'un agent de transfert et donc de sa constante de transfert

global Cs
Cs = 0.66; %Butanethiol

%%
%Déclaration des conditions initiales en [mol/m^2] rappel, le problème est
%infini dans 2 dimensions

global condInitM
condInitII  = 3e-2; %Amorçeur (arbitraire)
condInitM   = 9.4 ; %Monomère (fixée par les données de l'énoncé)

if(transfert == 1)
    condInitTrH = 5e-4; %Agent de transfert (arbitraire)
else
    condInitTrH = 0.0  ;
end

Ti = 295.65; % Condition initiale en température [K]

condInit    = [condInitII condInitM condInitTrH Ti]; %remise sous forme de vecteur pour la fonction ode15s

%%
% Partie Résolution du problème par ode15s
% Note, on écrit C pour la concentration, mais on résout également la
% température, qui sera contenu dans C(:,4)

n = 50000; %Nombre de pas de temps
tSpan  = linspace(0,86400,n); %Simulation sur 24 heures
[t, C] = ode15s(@(t,C)dCdt(t, C),tSpan,condInit); %Résolution des équations. On obtient des concentrations en [mol/m^2] et des températures en [K]
%ordre dans C = [amorçeur monomère agentTransfert température]

%%
% Partie calcul du Degré de polymérisation moyen instantané en nombre. On
% utilise les concentrations et températures obtenues par Ode15s pour
% calculer la probabilité de propagation : alpha

DPn = zeros(n,1);
ri = zeros(n,1);
for i=1:n
    
    x = 1 - (C(i,2)/condInitM); %Conversion
    T = C(i,4);
    
    if vitrification %Calcul de kp
        kp = calculKp(x,T);
    else
        kp = calculKp(0,T);
    end
        ks = Cs * kp; %Calcul de ks

    if trommsdorf % Calcul de kt
        kt = calculKt(x,T);
    else
        kt = calculKt(0,T);
    end

    if ricst %Calcul de la vitesse d'apparition des radicaux
        ri(i) = 8.36e-9; %Cas ri = cst
    else
        %Calcul de kII
        %AII = 4.0693e11 ; % Facteur pré-exponentielle
        %EII = 1.04192e5 ; % Energie d'activation
        %R   = 8.314     ;
        %kII = AII*exp(-EII/(R*T));% Calcul par la loi d'Arrhénius
        kII = 1.58e-7;
        ri(i) = 2*kII*C(1); % Cas ri = fct([I-I],T)
    end

    R = sqrt(ri(i)/kt); %Concentration totale en radicaux

    alpha = kp*R*C(i,2) / (kp*R*C(i,2) + ks*C(i,3)*R + kt*R^2); % Probabilité de propagation
    DPn(i) = 1/(1-alpha); % Dpn instantané moyen en nombre
    
end
%%
% Partie Graphes
if(0)
figure % Concentration en amorçeur [mol/m^2]
plot(t/3600,C(:,1));
xlabel('Temps [heure]');
xlim([0 24]);
ylabel('Concentration en amorçeur [mol/m^2]');

figure % Conversion en monomère [/]
plot(t/3600,1 - (C(:,2)/condInitM));
xlabel('Temps [heure]');
xlim([0 24]);
ylabel('Conversion en monomère [/]');

figure % Concentration en agent de tranfert [mol/m^2]
plot(t/3600, C(:,3));
xlabel('Temps [heure]');
xlim([0 24]);
ylabel('Concentration en agent de transfert [mol/m^2]');

figure % Température [°C]
plot(t/3600,C(:,4)-273.15);
xlabel('Temps [heure]');
xlim([0 24]);
ylabel('Température [°C]');

figure
plot(t/3600,DPn); % Degré de polymérisation instantanné moyen en nombre [/]
xlabel('Temps [heure]');
xlim([0 24]);
ylabel('Dpn [/]');
end

end


function Cd = dCdt (t, C)
global trommsdorf vitrification ricst temperature
global condInitM Cs 

x = 1 - (C(2)/condInitM);%Conversion du monomère
T = C(4);
%%
% Calcul des constantes cinétiques, en considérant certains effets ou non

if vitrification %Calcul de kp
    kp = calculKp(x,T);
else
    kp = calculKp(0,T);
end
ks = Cs * kp; %Calcul de ks

if trommsdorf % Calcul de kt
    kt = calculKt(x,T);
else
    kt = calculKt(0,T);
end

if ricst %Calcul de la vitesse d'apparition des radicaux
    ri = 8.36e-9; %Cas ri = cst
else
    %Calcul de kII
%     AII = 4.0693e11 ; % Facteur pré-exponentielle
%     EII = 1.04192e5 ; % Energie d'activation
%     R   = 8.314     ;
%     kII = AII*exp(-EII/(R*T));% Calcul par la loi d'Arrhénius
    kII = 1.58e-7;
    ri = 2*kII*C(1); % Cas ri = fct([I-I],T)
end

R = sqrt(ri/kt); %Concentration totale en radicaux

%%
% Calcul des vitesses de disparition des différents réactifs
% ordre dans C et Cd= [amorçeur monomère agentTransfert Température]
Cd(1) = -ri/2;
Cd(2) = -ri - kp*C(2)*R - ks*C(3)*R;
Cd(3) = -ks*C(3)*R;

%%
% Calcul de l'évolution de la température, si on ne considère pas
% d'impact, on fixe cette évolution à 0.0
if temperature
    cp = 114.1 + 6.83*T;
    DeltaH = -57800;
    rho = 0.94;
    Cd(4) = (DeltaH*Cd(2)) / (rho*cp);
else
    Cd(4) = 0.0;
end


Cd = Cd';


end

function kp = calculKp (x,T)
% Calcul de la constante cinétique associée à la propagation en fonction de
% la conversion et de la température. On approxime l'impact de la
% conversion par une partie constante jusqu'à x = 0.4 puis par un polynôme
% de degré 2 pour 0.4 < x < 1.0

Ep = 1.8283e4; %Energie d'ativation de la réaction de propagation

if(x < 0.4) % La conversion n'impact pas kp tant qu'elle ne dépasse pas 0.4
    Ap1 = 537412.19; % Facteur pré exponentiel pour cette partie, calculée sur base du graphe
    kp = Ap1 * exp(-Ep/(8.314*T)); %Calcul de kp
    
elseif(x >= 0.4)% Approximation de l'impact de la conversion par un polynome du 2ème degré
    Ap2 = 45165.57; % Facteur pré exponentiel pour cette partie, calculée sur base du graphe
    pol = [-10.7148 7.0289 1.4245]; %Coefficient du polynôme d'approximation
    kp = Ap2 * 10^( pol(1)*x^2 + pol(2)*x) * exp(-Ep/(8.314*T)); %Calcul de kp   
end
end


function kt = calculKt (x,T)
% Calcul de la constante cinétique associée à la terminaison en fonction de
% la conversion et de la température. On modèlise l'impact de la conversion
% par des polynomes du 1er degré pour 2 parties distinctes du graphe x<0.1
% et x>0.1

Et = 2.9442e3; % Energie d'activation de la réaction de terminaison

if(x<0.1)
    At1 = 1.47976e8; % Facteur pré exponentiel pour cette partie, calculée sur base du graphe
    kt = At1 * 10^(-1.818*x) * exp(-Et/(8.314*T)); % Calcul du kt
else
    At2 = 4.17054e8; % Facteur pré exponentiel pour cette partie, calculée sur base du graphe
    kt = At2 * 10^(-6.12*x)  * exp(-Et/(8.314*T)); % Calcul du kt
end

end

function[] = approximation()
%% 
% Cette fonction donne les coeffcients du polynôme de degré 2 ncessaires
% pour le calcul de la constante cinétique de propagation

y = [2.5 2.3 1.8 1.2 0.0];
conversion = [0.4  0.5  0.6  0.68 0.82];
pol = polyfit(conversion,y,2);

end


