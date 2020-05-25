%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%             TP1 : Traitement du signal            %%%%%%%%%%%%%
%%%%%%%%%%%             Auteur : THUET Quentin                %%%%%%%%%%%%%
%%%%%%%%%%%                      CAMARA Ababacar              %%%%%%%%%%%%%
%%%%%%%%%%%             Filiere: 1Sn.                         %%%%%%%%%%%%%
%%%%%%%%%%%             Groupe F.                             %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear;
close all;
taille_ecran = get(0,'ScreenSize');
L = taille_ecran(3);
H = taille_ecran(4);



%% 2.2 Premiere Chaine : " Chaine de reference "



%% 2.2.1

% Duree symbole en nombre d'echantillons (Ts = Ns*Te)
Ns = 10;
% Nombre de bits generes
nb_bits = 10100;
% Generation des bits
bits = randi([0 1],1,nb_bits);
% Mapping binaire a moyenne nulle : 0 ->-1, 1->1
Symboles = 2 * bits - 1;
% Generation de la suite de diracs ponderee par les symboles (surechantillonnage)
Suite_diracs = kron(Symboles, [1 zeros(1,Ns-1)]);
% Generation de la reponse impulsionnelle du filtre de mise en forme (NRZ)
h = ones(1,Ns);
% Filtrage de mise en forme
x = filter(h,1,Suite_diracs);
% Affichage du signal genere
figure(1);clf
plot(x);
title('Signal emis')
xlabel('bit')
ylabel('amplitude')
axis([0 200 -1.5 1.5]);
% Calcul de la DSP du signal par periodogramme
DSP_x = (1/length(x))*abs(fft(x,2^nextpow2(length(x)))).^2;
% Affichage de la DSP du signal genere
figure(2);clf
plot(linspace(-1,1,length(DSP_x)), fftshift(DSP_x));
title('DSP du signal emis')
xlabel('frequences normalisees')
ylabel('amplitude')

%% 2.2.2.

%% 2.2.2 - a génération du signal en sortie du filtre de reception
% calcul de g(t)
hr = ones(1,Ns);
g  = conv(h,hr);

figure(3); clf
plot(g);
axis([0 Ns*2 0 25]);
title('g = h * h_r')
xlabel('t (secondes)')
ylabel('g(t)')

%% 2.2.2 - b  diagramme de l'oeil

% Generation du signal recu
z = conv(x,hr);
figure(4);clf
plot(z);
axis([0 200 -(Ns+5) (Ns+5)]);
title('Signal en sortie de filtre de reception')
xlabel('bit')
ylabel('amplitude')

% Diagramme de l'oeil
eyediagram(z(Ns*2:end-(Ns+1)),Ns*2,Ns*2);

%% 2.2.2 - c Prise de decision et calcul du taux d'erreur binaire

% Sous-echantillonage (instant optimaux)
t0 = Ns;
z_ech = z(t0:Ns:end);
% Symboles decides
symboles_decides = sign(z_ech);
% Demapping
bits_decides = (1/2) * (symboles_decides + 1);
nb_erreurs = length(find(bits ~= bits_decides));
teb = nb_erreurs/length(bits);
disp(['Le TEB pour la chaine 1 non bruitee est : ' num2str(teb)]);


%% 2.2.3  Implantation du canal bruite

z_bruit = bruitage(x,hr,Ns,2,2);
% Sous-echantillonage (instant optimaux)
z_ech_bruit = z_bruit(Ns:Ns:end);

figure(6);clf
plot(z_bruit);
axis([0 200 -(Ns+5) (Ns+5)]);
title('Signal en sortie de filtre de reception bruite')
xlabel('bit')
ylabel('amplitude')

% Diagramme de l'oeil
eyediagram(z_bruit(Ns*2:end-(Ns+1)),Ns*2,Ns*2);


%% 2.2.4 Superposition du TEB simule et du TEB theorique

tab_En_b_N0 = [0 1 2 3 4 5 6 ];
tab_teb_theorique = qfunc(sqrt(2*10.^(tab_En_b_N0/10)));

tab_teb_simule=zeros(1,7);
for i = 0:6
    
    z_bruit = bruitage(x,hr,Ns,2,i);
    % Sous-echantillonage (instant optimaux)
    z_ech_bruit = z_bruit(Ns:Ns:end);
    bits_decides_bruit = decision2aire(z_ech_bruit);
    teb = tauxErreur(bits,bits_decides_bruit);
    
    tab_teb_simule(i+1)=teb;
end

figure(8); clf
semilogy(tab_En_b_N0,tab_teb_theorique,tab_En_b_N0,tab_teb_simule);
title('Comparaison TEB simule et TEB theorique')
xlabel('RSB (dB)')
ylabel('TEB')
lg = legend('TEB theorique', ...
    'TEB simule', ...
    'Location','Best');
set(lg,'Interpreter','Latex');



%% 3.2 Deuxieme chaine : " Impact du choix du filtre de reception "


%% 3.2.1

%% 3.2.1.a Trace du signal en sortie du filtre de reception

% Filtre de reception
hr2 = ones(1,Ns);
hr2(Ns/2+1:end) = 0;
% Generation du signal en sortie du filtre de reception
z2 = conv(x,hr2);
figure(9);clf
plot(z2);
axis([0 200 -(Ns+5) (Ns+5)]);
title('Signal en sortie de filtre de reception')
xlabel('bit')
ylabel('amplitude')

%% 3.2.1.b Trace du diagramme de l'oeil

eyediagram(z2(Ns*2:end-(Ns+1)),Ns*2,Ns*2);

%% 3.2.1.c Calcul du TEB

t0 = Ns;
z_ech = z2(t0:Ns:end);
% Symboles decides
symboles_decides = sign(z_ech);
% Demapping
bits_decides = (1/2) *( symboles_decides + 1);
nb_erreurs = length(find(bits ~= bits_decides));
teb2 = nb_erreurs/length(bits);
disp(['Le TEB pour la chaine 2 non bruitee est : ' num2str(teb2)]);

%% 3.2.2 Implantation de la chaine avec bruit

z2_bruit = bruitage(x,hr2,Ns,2,15);

figure(11);clf
plot(z2_bruit);
axis([0 200 -(Ns+5) (Ns+5)]);
title('Signal en sortie de filtre de reception bruite')
xlabel('bit')
ylabel('amplitude')
eyediagram(z2_bruit((Ns*2):end-(Ns+1)),Ns*2,Ns*2);

%% 3.3.3  Trace des TEB simule et theorique

tab_En_b_N0 = [0 1 2 3 4 5 6 ];
tab_teb_theorique2 = qfunc(sqrt(10.^(tab_En_b_N0/10)));

tab_teb_simule2 = zeros(1,7);

for i = 0:6
    z2_bruit = bruitage(x,hr2,Ns,2,i);
    % Sous-echantillonage (instant optimaux)
    z2_ech_bruit = z2_bruit(Ns:Ns:end);
    bits_decides_bruit = decision2aire(z2_ech_bruit);
    teb2 = tauxErreur(bits,bits_decides_bruit);
    
    tab_teb_simule2(i+1) = teb2;
end

figure(13); clf
semilogy(tab_En_b_N0,tab_teb_theorique2,tab_En_b_N0,tab_teb_simule2);
title('Comparaison TEB simule et TEB theorique')
xlabel('RSB (dB)')
ylabel('TEB')
lg = legend('TEB theorique', ...
    'TEB simule', ...
    'Location','Best');
set(lg,'Interpreter','Latex');

%% 3.3.4  Comparaison avec la chaine de reference

figure(14); clf
semilogy(tab_En_b_N0,tab_teb_simule,tab_En_b_N0,tab_teb_simule2);
title('Comparaison TEB de la premiere chaine et TEB de la deuxieme chaine')
xlabel('RSB (dB)')
ylabel('TEB')
lg = legend('TEB simule chaine reference', ...
    'TEB simule seconde chaine', ...
    'Location','Best');
set(lg,'Interpreter','Latex');

%% 3.3.5 Comparaison de l'efficacite spectrale des deux chaines de transmission

% Calcul de la DSP du signal par periodogramme pour le signal de reference
DSP_x = (1/length(x))*abs(fft(x,2^nextpow2(length(x)))).^2;
% Calcul de la DSP du signal etudie
DSP_x_etudie = (1/length(x))*abs(fft(x,2^nextpow2(length(x)))).^2;

figure(15);clf
plot(linspace(-1,1,length(DSP_x)), fftshift(DSP_x),linspace(-1,1,length(DSP_x_etudie)), fftshift(DSP_x_etudie));
title('Comparaison DSP de la premiere chaine et DSP de la deuxieme chaine')
xlabel('frequences normalisees')
ylabel('amplitude')
lg = legend('DSP chaine de reference', ...
    'DSP seconde chaine', ...
    'Location','Best');
set(lg,'Interpreter','Latex');



%% 4.2 Troisieme chaine : Impact du choix du filtre de mise en forme et du canal de propagation



%% 4.2.1 Initialisation

Fe = 12000;
Rs = 3000;
alpha = 0.5;    % Roll-off factor
Ns = Fe / Rs;
% Generation de la suite de Diracs ponderes par les symboles (surechantillonnage)
Suite_diracs2 = kron(Symboles, [1 zeros(1,Ns-1)]);

%% 4.2.2

%% 4.2.2-a Tracer le signal en sortie du filtre de reception

% Filtre de mise en forme (cosinus sureleve)
Nsym = Ns;           % Filter span in symbol durations
sampsPerSym = Ns;    % Upsampling factor
h2 = rcosdesign(alpha, Ns, sampsPerSym, 'sqrt');
% Signal en sortie du filtre de mise en forme
x2 = conv(Suite_diracs2,h2);
figure(16);clf
plot(x2);
axis([0 200 -(Ns+5) (Ns+5)]);
title('Signal en sortie du filtre de mise en forme');
xlabel('bit');
ylabel ('amplitude');
% Filtre de reception
hr3 = rcosdesign(alpha, Ns, sampsPerSym, 'sqrt');
% Signal en sortie du filtre de reception
z3 = conv (x2 ,hr3);
figure(17);clf
plot(z3);
axis([0 200 -(Ns+5) (Ns+5)]);
title('Signal en sortie du filtre de reception bruite')
xlabel('bit')
ylabel('amplitude')

%% 4.2.2-b

% Diagramme de l'oeil
eyediagram(z3(Ns+1:end-Ns),2*Ns,Ns*2);

%% 4.2.2-c Prise de decision

% Sous-echantillonage (instants optimaux)
t0 = length(h2);
z_ech = z3(t0:Ns:end-t0);
% Symboles decides
symboles_decides = sign(z_ech);
% Demapping
bits_decides = (1/2) *(symboles_decides + 1);
nb_erreurs = length(find(bits ~= bits_decides));
teb = nb_erreurs/length(bits);
disp(['Le TEB pour la chaine 3 non bruitee est : ' num2str(teb)]);

%% 4.2.3 Implantation du bruit

z_bruit = bruitage(x2,h2,Ns,2,1);

%% 4.2.4 Comparaison TEB simule et TEB theorique de la chaine de transmission avec le racine de cosinus sureleve

tab_En_b_N0 = [0 1 2 3 4 5 6 ];
tab_teb_theorique3 = qfunc(sqrt(2*10.^(tab_En_b_N0/10)));

tab_teb_simule3 = zeros(1,7);
for i = 0:6
    z_bruit = bruitage(x2,h2,Ns,2,i);
    % Sous-echantillonage (instant optimaux)
    z_ech_bruit = z_bruit(t0:Ns:end-t0);
    bits_decides_bruit = decision2aire(z_ech_bruit);
    teb3=tauxErreur(bits,bits_decides_bruit);
    
    tab_teb_simule3(i+1) = teb3;
end

figure(19); clf
semilogy(tab_En_b_N0,tab_teb_theorique3,tab_En_b_N0,tab_teb_simule3);
title('Comparaison TEB simule et TEB theorique')
xlabel('RSB (dB)')
ylabel('TEB')
lg = legend('TEB theorique', ...
    'TEB simule', ...
    'Location','Best');
set(lg,'Interpreter','Latex');

%% 4.2.5 Comparaison du taux d'erreur binaire obtenu avec la racine de cosinus sureleve avec celui de la chaine de reference

figure(20); clf
semilogy(tab_En_b_N0,tab_teb_simule,tab_En_b_N0,tab_teb_simule3);
title('Comparaison TEB premiere chaine et TEB troisieme chaine')
xlabel('RSB (dB)')
ylabel('TEB')
lg = legend('TEB simule chaine reference', ...
    'TEB simule raised cosinus filter', ...
    'Location','Best');
set(lg,'Interpreter','Latex');

%% 4.2.6 Comparaison de l'efficacite spectrale

% Calcul de la DSP du signal par periodogramme pour le signal de reference
DSP_x = (1/length(x))*abs(fft(x,2^nextpow2(length(x)))).^2;
% Calcul de la DSP du signal etudie
DSP_x_etudie = (1/length(x2))*abs(fft(x2,2^nextpow2(length(x2)))).^2;

figure(21);clf
plot(linspace(-1,1,length(DSP_x)), fftshift(DSP_x),linspace(-1,1,length(DSP_x_etudie)), fftshift(DSP_x_etudie));

title('Comparaison DSP de la premiere chaine et DSP de la troisieme chaine')
xlabel('frequences normalisees')
ylabel('amplitude')
lg = legend('DSP reference', ...
    'DSP racine cosinus sureleve.', ...
    'Location','Best');
set(lg,'Interpreter','Latex');

%% 4.2.7 Chaine de transmission avec un passage dans un canal de bande BW

%% 4.2.7-a BW = 1500 Hz

N = 60;
fc1 = 1500;
morceau = (-N:N)/Fe;
h_pb = 2*fc1/Fe*sinc(2*fc1*morceau);
x2reduit = conv (h_pb,x2);

z4 = conv (x2reduit ,hr3);
% Diagramme de l'oeil
eyediagram(z4(Ns+1:end-Ns),2*Ns,Ns*2);

%% 4.2.7-b BW = 3000 Hz

fc2 = 3000;
morceau = (-N:N)/Fe;
h_pb = 2*fc2/Fe*sinc(2*fc2*morceau);
x2reduit2 = conv (h_pb,x2);
z5 = conv (x2reduit2 ,hr3);
% Diagramme de l'oeil
eyediagram(z5(Ns+1:end-Ns),2*Ns,Ns*2);



%% 5.2 Quatrieme Chaine : " Impact du choix du mapping "



%% 5.2.1

% Mapping
Ns = 10*2;
symboles2 = (2*bi2de(reshape(bits,2,length(bits)/2).')-3).';
Suite_diracs3 = kron(symboles2, [1 zeros(1,Ns-1)]);
h = ones(1,Ns);

%% 5.2.1-a

% Signal en sortie du filtre d'emission
x5 = conv(Suite_diracs3,h);
figure(24);clf;
plot(x5);
axis([0 400 -(5) (5)]);
title('Signal emis')
xlabel('bit')
ylabel('amplitude')
% Calcul de la DSP du signal par periodogramme
DSP_x = (1/length(x5))*abs(fft(x5,2^nextpow2(length(x5)))).^2;
% Affichage de la DSP du signal genere
figure(25);clf
plot(linspace(-1,1,length(DSP_x)), fftshift(DSP_x));
title('DSP du signal emis')
xlabel('frequences normalisees')
ylabel('amplitude')

%% 5.2.1-b Comparaison de l'efficacite spectrale

% Calcul de la DSP du signal par periodogramme pour le signal de reference
DSP_x = (1/length(x))*abs(fft(x,2^nextpow2(length(x)))).^2/80;
% Calcul de la DSP du signal etudie
DSP_x_etudie = (1/length(x5))*abs(fft(x5,2^nextpow2(length(x5)))).^2/600;
figure(26);clf
plot(linspace(-1,1,length(DSP_x)), fftshift(DSP_x),linspace(-1,1,length(DSP_x_etudie)), fftshift(DSP_x_etudie));
title('Comparaison DSP de la premiere chaine et DSP de la quatrieme chaine')
xlabel('frequences normalisees')
ylabel('amplitude')
lg = legend('DSP reference', ...
    'DSP racine cosinus sureleve.', ...
    'Location','Best');
set(lg,'Interpreter','Latex');

%% 5.2.1-c Tracer le diagramme de l'oeil

% Diagramme de l'oeil
z6 = conv(x5,h);
eyediagram(z6(Ns*2:end-Ns),2*Ns,Ns*2);

%% 5.2.1-d Prise de decision

% Sous-echantillonage (instant optimaux)
t0 = Ns;
z_ech = z6(t0:Ns:end-Ns);
% Symboles decides
decisionMoins3 = (z_ech <= (-2*Ns))*(-3);
decisionMoins1 = ((-2*Ns<z_ech)&(z_ech<=0))*(-1);
decisionPlus1 = (0<z_ech)&(z_ech<=2*Ns);
decisionPlus3 = (z_ech >= (2*Ns))*(3);
symboles_decides= decisionMoins3+decisionMoins1+decisionPlus1+decisionPlus3;
% Demapping
bits_decides = reshape(de2bi((symboles_decides+3)/2).',1,length(bits));
nb_erreurs = length(find(bits ~= bits_decides));
teb = nb_erreurs/length(bits);
disp(['Le TEB pour la chaine 4 non bruitee est : ' num2str(teb)]);

%% 5.2.2 Taux d'erreur symbole simule en fonction du rapport signal sur bruit

tab_En_b_N0 = [0 1 2 3 4 5 6 ];
tab_tes_simule5 = zeros(1,7);
for i = 0:6
    z_bruit = bruitage(x5,h,Ns,4,i);
    % Sous-echantillonnage (instants optimaux)
    z_ech_bruit = z_bruit(t0:Ns:end-Ns);
    % Symboles decides
    symboles_decides= decision4aire(z_ech_bruit,Ns);
    % Demapping
    tes5=tauxErreur(symboles2,symboles_decides);
    
    tab_tes_simule5(i+1) = tes5;
end

figure(28); clf
semilogy(tab_En_b_N0,tab_tes_simule5);
title('TES')
xlabel('RSB (dB)')
ylabel('TES')

%% 5.2.3 Comparaison TES simule et TES theorique

figure(29); clf
tab_tes_theorique5 = (3/2)* qfunc(sqrt((4/5)*10.^(tab_En_b_N0/10)));
semilogy(tab_En_b_N0,tab_tes_theorique5,tab_En_b_N0,tab_tes_simule5);
title('Comparaison TES simule et TES theorique')
xlabel('RSB (dB)')
ylabel('TES')
lg = legend('TES simule chaine theorique', ...
    'TES simule', ...
    'Location','Best');
set(lg,'Interpreter','Latex');

%% 5.2.4 Trace du taux d'erreur binaire

figure(30); clf
tab_teb_simule5 = zeros(1,7);
for i = 0:6
    z_bruit = bruitage(x5,h,Ns,4,i);
    % Sous-echantillonage (instant optimaux)
    z_ech_bruit = z_bruit(t0:Ns:end-Ns);
    % Symboles decides
    symboles_decides = decision4aire(z_ech_bruit,Ns);
    % Demapping
    bits_decides = reshape(de2bi((symboles_decides+3)/2).',1,length(bits));
    teb5=tauxErreur(bits,bits_decides);
    
    tab_teb_simule5(i+1) = teb5;
end
semilogy(tab_En_b_N0,tab_teb_simule5);

%% 5.2.4 Comparaison du TEB theorique au TEB simule

tab_En_b_N0 = [0 1 2 3 4 5 6 ];
tab_teb_theorique5 = (1/2)*tab_tes_theorique5;
semilogy(tab_En_b_N0,tab_teb_theorique5,tab_En_b_N0,tab_teb_simule5);
title('Comparaison TEB simule et TEB theorique')
xlabel('RSB (dB)')
ylabel('TEB')
lg = legend('TEB  theorique (Mapping de Gray)', ...
    'TEB simule (Mapping Naturel)', ...
    'Location','Best');
set(lg,'Interpreter','Latex');

%% Fin du TP. Merci de la lecture.





