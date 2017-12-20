function  [DPn,C,ri] = Cinetique2017(a,b,c,d,e)


%%
%Cette section contient toutes les options de simulations choisies, 0 pour non, 1
%pour oui. On a, dans l'ordre, l'effet Trommsdorf, l'effet de
%vitrification, l'ajout d'un agent de transfert, l'�puisement de l'amor�eur
%et l'impact de la temp�rature. Ce variables sont utilis�es plus tard dans
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
%D�claration des conditions initiales en [mol/m^2] rappel, le probl�me est
%infini dans 2 dimensions

global condInitM
condInitII  = 3e-2; %Amor�eur (arbitraire)
condInitM   = 9.4 ; %Monom�re (fix�e par les donn�es de l'�nonc�)

if(transfert == 1)
    condInitTrH = 5e-4; %Agent de transfert (arbitraire)
else
    condInitTrH = 0.0  ;
end

Ti = 295.65; % Condition initiale en temp�rature [K]

condInit    = [condInitII condInitM condInitTrH Ti]; %remise sous forme de vecteur pour la fonction ode15s

%%
% Partie R�solution du probl�me par ode15s
% Note, on �crit C pour la concentration, mais on r�sout �galement la
% temp�rature, qui sera contenu dans C(:,4)

n = 50000; %Nombre de pas de temps
tSpan  = linspace(0,86400,n); %Simulation sur 24 heures
[t, C] = ode15s(@(t,C)dCdt(t, C),tSpan,condInit); %R�solution des �quations. On obtient des concentrations en [mol/m^2] et des temp�ratures en [K]
%ordre dans C = [amor�eur monom�re agentTransfert temp�rature]

%%
% Partie calcul du Degr� de polym�risation moyen instantan� en nombre. On
% utilise les concentrations et temp�ratures obtenues par Ode15s pour
% calculer la probabilit� de propagation : alpha

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
        %AII = 4.0693e11 ; % Facteur pr�-exponentielle
        %EII = 1.04192e5 ; % Energie d'activation
        %R   = 8.314     ;
        %kII = AII*exp(-EII/(R*T));% Calcul par la loi d'Arrh�nius
        kII = 1.58e-7;
        ri(i) = 2*kII*C(1); % Cas ri = fct([I-I],T)
    end

    R = sqrt(ri(i)/kt); %Concentration totale en radicaux

    alpha = kp*R*C(i,2) / (kp*R*C(i,2) + ks*C(i,3)*R + kt*R^2); % Probabilit� de propagation
    DPn(i) = 1/(1-alpha); % Dpn instantan� moyen en nombre
    
end
%%
% Partie Graphes
if(0)
figure % Concentration en amor�eur [mol/m^2]
plot(t/3600,C(:,1));
xlabel('Temps [heure]');
xlim([0 24]);
ylabel('Concentration en amor�eur [mol/m^2]');

figure % Conversion en monom�re [/]
plot(t/3600,1 - (C(:,2)/condInitM));
xlabel('Temps [heure]');
xlim([0 24]);
ylabel('Conversion en monom�re [/]');

figure % Concentration en agent de tranfert [mol/m^2]
plot(t/3600, C(:,3));
xlabel('Temps [heure]');
xlim([0 24]);
ylabel('Concentration en agent de transfert [mol/m^2]');

figure % Temp�rature [�C]
plot(t/3600,C(:,4)-273.15);
xlabel('Temps [heure]');
xlim([0 24]);
ylabel('Temp�rature [�C]');

figure
plot(t/3600,DPn); % Degr� de polym�risation instantann� moyen en nombre [/]
xlabel('Temps [heure]');
xlim([0 24]);
ylabel('Dpn [/]');
end

end


function Cd = dCdt (t, C)
global trommsdorf vitrification ricst temperature
global condInitM Cs 

x = 1 - (C(2)/condInitM);%Conversion du monom�re
T = C(4);
%%
% Calcul des constantes cin�tiques, en consid�rant certains effets ou non

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
%     AII = 4.0693e11 ; % Facteur pr�-exponentielle
%     EII = 1.04192e5 ; % Energie d'activation
%     R   = 8.314     ;
%     kII = AII*exp(-EII/(R*T));% Calcul par la loi d'Arrh�nius
    kII = 1.58e-7;
    ri = 2*kII*C(1); % Cas ri = fct([I-I],T)
end

R = sqrt(ri/kt); %Concentration totale en radicaux

%%
% Calcul des vitesses de disparition des diff�rents r�actifs
% ordre dans C et Cd= [amor�eur monom�re agentTransfert Temp�rature]
Cd(1) = -ri/2;
Cd(2) = -ri - kp*C(2)*R - ks*C(3)*R;
Cd(3) = -ks*C(3)*R;

%%
% Calcul de l'�volution de la temp�rature, si on ne consid�re pas
% d'impact, on fixe cette �volution � 0.0
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
% Calcul de la constante cin�tique associ�e � la propagation en fonction de
% la conversion et de la temp�rature. On approxime l'impact de la
% conversion par une partie constante jusqu'� x = 0.4 puis par un polyn�me
% de degr� 2 pour 0.4 < x < 1.0

Ep = 1.8283e4; %Energie d'ativation de la r�action de propagation

if(x < 0.4) % La conversion n'impact pas kp tant qu'elle ne d�passe pas 0.4
    Ap1 = 537412.19; % Facteur pr� exponentiel pour cette partie, calcul�e sur base du graphe
    kp = Ap1 * exp(-Ep/(8.314*T)); %Calcul de kp
    
elseif(x >= 0.4)% Approximation de l'impact de la conversion par un polynome du 2�me degr�
    Ap2 = 45165.57; % Facteur pr� exponentiel pour cette partie, calcul�e sur base du graphe
    pol = [-10.7148 7.0289 1.4245]; %Coefficient du polyn�me d'approximation
    kp = Ap2 * 10^( pol(1)*x^2 + pol(2)*x) * exp(-Ep/(8.314*T)); %Calcul de kp   
end
end


function kt = calculKt (x,T)
% Calcul de la constante cin�tique associ�e � la terminaison en fonction de
% la conversion et de la temp�rature. On mod�lise l'impact de la conversion
% par des polynomes du 1er degr� pour 2 parties distinctes du graphe x<0.1
% et x>0.1

Et = 2.9442e3; % Energie d'activation de la r�action de terminaison

if(x<0.1)
    At1 = 1.47976e8; % Facteur pr� exponentiel pour cette partie, calcul�e sur base du graphe
    kt = At1 * 10^(-1.818*x) * exp(-Et/(8.314*T)); % Calcul du kt
else
    At2 = 4.17054e8; % Facteur pr� exponentiel pour cette partie, calcul�e sur base du graphe
    kt = At2 * 10^(-6.12*x)  * exp(-Et/(8.314*T)); % Calcul du kt
end

end

function[] = approximation()
%% 
% Cette fonction donne les coeffcients du polyn�me de degr� 2 ncessaires
% pour le calcul de la constante cin�tique de propagation

y = [2.5 2.3 1.8 1.2 0.0];
conversion = [0.4  0.5  0.6  0.68 0.82];
pol = polyfit(conversion,y,2);

end


