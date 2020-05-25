clear
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3.2 Utilisation de la chaine passe-bas equivalente pour le calcul et l'estimation du TEB %%
%%                                                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3.2.1 Implantation de la chaine sur frequence porteuse %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alpha = 0.35;       % Roll-off factor
span = 20;          % Roll-off span factor
fp = 2000;          % Frequence porteuse
fe = 10000;         % Frequence d'echantillonage
Te = 1 / fe;        % Periode d'echantillonage
Rs = 1000;          % Debit symbole
Ts = 1 / Rs;        % Duree symbole
Ns = fe / Rs;       % Facteur de surechantillonage
nb_bits = 10100;    % Taille de l'information a transmettre 
N = 100;            % Ordre du filtre passe-bas

%% 3.2.1.1 Trace du signal genere et du signal transmis

% Generation de l'information a transmettre
bits = randi([0 1],1,nb_bits);
symboles = 2*bits - 1;

% Mapping
a = symboles(1:2:end);
b = symboles(2:2:end);

% Surechantillonage
a_temp = kron(a, [1 zeros(1,Ns-1)]);
b_temp = kron(b, [1 zeros(1,Ns-1)]);

% Introduction de zeros en debut et fin de chaine pour anticiper le retard
% du au filtre de mise en forme
a = zeros(1, length(a_temp)+span*Ns);
b = zeros(1, length(b_temp)+span*Ns);
a(span/2*Ns+1:end-span/2*Ns) = a_temp;
b(span/2*Ns+1:end-span/2*Ns) = b_temp;

% Filtrage de mise en forme 
h = rcosdesign(alpha, span, Ns, 'sqrt');

I_retard = conv(a,h);
Q_retard = conv(b,h);

% Elimination du retard
I = I_retard(span*Ns+1:end-span*Ns);
Q = Q_retard(span*Ns+1:end-span*Ns);

t = linspace(0, nb_bits/2*Ts, nb_bits/2*Ns); 

figure(1); clf
subplot(211);
plot(t,I);
title('Signal genere sur la voie en phase');
xlabel('t (s)');
ylabel('I(t)');
subplot(212);
plot(t,Q);
title('Signal genere sur la voie en quadrature');
xlabel('t (s)');
ylabel('Q(t)');

% Creation de l'enveloppe complexe
xe = I + 1j*Q;

% Transposition de frequence
trans = exp(1j*2*pi*fp*t);

x = real(xe.*trans);

figure(2); clf
plot(t,x);
title('Signal transmis sur frequence porteuse');
xlabel('t (s)');
ylabel('x(t)');

%% 3.2.1.2 Trace de la DSP du signal transmis

DSP_x = (1/length(x))*abs(fft(x,...
    2^nextpow2(length(x)))).^2;

f = linspace(-fe/2,fe/2,length(DSP_x));

figure(3); clf
plot(f,fftshift(DSP_x));
title('DSP du signal transmis')
xlabel('frequence (Hz)')
ylabel('amplitude')

%% 3.2.1.3 Chaine complete sans bruit

% Retour en bande de base
zb = bande_de_base(x, t, fp, fe, N);

% Filtrage de reception
z = filtrage_reception(zb, span, Ns, h);

% Echantillonage, decision et demapping
bits_decides = decision(z, Ns, nb_bits, 'psk', 4, false);

disp(['TEB chaine sur frequence porteuse sans bruit = ' num2str(TEB(bits_decides, bits, nb_bits)) '%']);

%% 3.2.1.4 Ajout du bruit

RSB = 0:6;
teb_bruit = zeros(1,length(RSB));

for i = RSB
    zb = bande_de_base(bruitage(x, Ns, 4, i, false), t, fp, fe, N);
    z = filtrage_reception(zb, span, Ns, h);
    bits_decides = decision(z, Ns, nb_bits, 'psk', 4, false);
    teb_bruit(i+1) = TEB(bits_decides, bits, nb_bits);
end

figure(4); clf
plot(RSB, teb_bruit);
title('TEB pour differentes valeurs de bruit')
xlabel('RSB')
ylabel('TEB')

%% 3.2.1.5 Comparaison TEB théorique et TEB simule

figure(5); clf
plot(RSB, qfunc(sqrt(2*10.^(RSB/10))), RSB, teb_bruit);
title('TEB theorique et simule pour differentes valeurs de bruit')
xlabel('RSB')
ylabel('TEB')
lg = legend('TEB theorique', ...
    'TEB simule', ...
    'Location','Best');
set(lg,'Interpreter','Latex');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3.2.2 Implantation de la chaine passe-bas equivalente%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


teb_bruit_porteuse = teb_bruit;

%% 3.2.1.1 Trace du signal genere et du signal transmis

figure(6); clf
subplot(211);
plot(t,I);
title('Signal genere sur la voie en phase');
xlabel('t (s)');
ylabel('I(t)');
subplot(212);
plot(t,Q);
title('Signal genere sur la voie en quadrature');
xlabel('t (s)');
ylabel('Q(t)');

%% 3.2.1.2 Trace de la DSP du signal transmis

DSP_xe = (1/length(xe))*abs(fft(xe,...
    2^nextpow2(length(xe)))).^2;

f = linspace(-fe/2,fe/2,length(DSP_xe));

figure(7); clf
plot(f,fftshift(DSP_xe));
title('DSP du signal transmis')
xlabel('frequence (Hz)')
ylabel('amplitude')

%% 3.2.2.3 Chaine complete sans bruit

% Filtrage de reception
z = filtrage_reception(xe, span, Ns, h);

% Echantillonage, decision et demapping
bits_decides = decision(z, Ns, nb_bits, 'psk', 4, false);

disp(['TEB chaine equivalente sans bruit = ' num2str(TEB(bits_decides, bits, nb_bits)) '%']);

%% 3.2.2.4 Ajout du bruit

RSB = 0:6;
teb_bruit = zeros(1,length(RSB));

for i = RSB
    z_bruit = filtrage_reception(bruitage(xe, Ns, 4, i, true), span, Ns, h);
    bits_decides = decision(z_bruit, Ns, nb_bits, 'psk', 4, false);
    teb_bruit(i+1) = TEB(bits_decides, bits, nb_bits);
end

figure(8); clf
plot(RSB, teb_bruit);
title('TEB pour differentes valeurs de bruit')
xlabel('RSB')
ylabel('TEB')

%% 3.2.2.5 Trace des constellations en sortie de mapping et de l'echantillonneur

scatterplot(a_temp(1:Ns:end) + 1j * b_temp(1:Ns:end));
hold on
title('xe');
hold off
scatterplot(z(1:Ns:end));
hold on
title('z');
hold off

%% 3.2.2.6 Comparaison TEB simule de la chaine sur freq porteuse et de la chaine equivalente

figure(11); clf
plot(RSB, teb_bruit_porteuse, RSB, teb_bruit);
title('TEB simule de la chaine sur freq porteuse et de la chaine equivalente')
xlabel('RSB')
ylabel('TEB')
lg = legend('TEB freq proteuse', ...
    'TEB chaine equivalente', ...
    'Location','Best');
set(lg,'Interpreter','Latex');


