clear
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4.2 Comparaison de modulations sur frequence porteuse %%
%%                                              %%%%%%%%%%%
%% 4.2.1 Etude de chaque chaine de transmission %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alpha = 0.5;        % Roll-off factor
span = 20;          % Roll-off span factor
fe = 48000;         % Frequence d'echantillonage
Te = 1 / fe;        % Periode d'echantillonage
Rb = 48000;         % Debit binaire
M = [4, 4, 8, 16];  % Nombre de bits par symbole
Rs = Rb ./ log2(M); % Debit symbole
Ts = 1 ./ Rs;       % Duree symbole
Ns = fe ./ Rs;      % Facteur de surechantillonage
nb_bits = 12000;    % Taille de l'information a transmettre 

% Generation de l'information a transmettre
bits = randi([0 1],1,nb_bits);

%% 4.2.1.1 Chaines completes sans bruit

% Information dans les nouvelles bases
bits_M4 = 2*bits(1:2:end)+bits(2:2:end);
bits_M8 = 4*bits(1:3:end)+2*bits(2:3:end)+bits(3:3:end);
bits_M16 = 8*bits(1:4:end)+4*bits(2:4:end)+2*bits(3:4:end)+bits(4:4:end);

% Calcul des symboles
symboles_4ask = bits_M4*2-3;
symboles_4ask = symboles_4ask + 2*(symboles_4ask == 1) - 2*(symboles_4ask == 3); % Mapping de gray
symboles_qpsk = pskmod(bits_M4,4,pi/4,'gray');
symboles_8psk = pskmod(bits_M8, 8,pi/4,'gray');
symboles_16qam = qammod(bits_M16,16,'gray'); 

% Surechantillonage, mise en forme et reception
[z_4ask, xe_4ask, h_4ask] = emission_reception(symboles_4ask, Ns(1), alpha, span);
[z_qpsk, xe_qpsk, h_qpsk] = emission_reception(symboles_qpsk, Ns(2), alpha, span);
[z_8psk, xe_8psk, h_8psk] = emission_reception(symboles_8psk, Ns(3), alpha, span);
[z_16qam, xe_16qam, h_16qam] = emission_reception(symboles_16qam, Ns(4), alpha, span);

% Echantillonage, decision et demapping
bits_decides_4ask = decision(z_4ask, Ns(1), nb_bits, 'ask', 4, true);
bits_decides_qpsk = decision(z_qpsk, Ns(2), nb_bits, 'psk', 4, true);
bits_decides_8psk = decision(z_8psk, Ns(3), nb_bits, 'psk', 8, true);
bits_decides_16qam = decision(z_16qam, Ns(4), nb_bits, 'qam', 16, true);   

disp(['TEB chaine 4-ASK sans bruit = ' num2str(TEB(bits_decides_4ask, bits, nb_bits)) '%']);
disp(['TEB chaine QPSK sans bruit = ' num2str(TEB(bits_decides_qpsk, bits, nb_bits)) '%']);
disp(['TEB chaine 8-PSK sans bruit = ' num2str(TEB(bits_decides_8psk, bits, nb_bits)) '%']);
disp(['TEB chaine 16-QAM sans bruit = ' num2str(TEB(bits_decides_16qam, bits, nb_bits)) '%']);

%% 4.2.1.2 Ajout du bruit

RSB = 0:6;
teb_bruit = zeros(4,length(RSB));

for i = RSB
    z_4ask = filtrage_reception(bruitage(xe_4ask, Ns(1), 4, i, false), span, Ns(1), h_4ask);
    z_qpsk = filtrage_reception(bruitage(xe_qpsk, Ns(2), 4, i, true), span, Ns(2), h_qpsk);
    z_8psk = filtrage_reception(bruitage(xe_8psk, Ns(3), 8, i, true), span, Ns(3), h_8psk);
    z_16qam = filtrage_reception(bruitage(xe_16qam, Ns(4), 16, i, true), span, Ns(4), h_16qam);
    bits_decides_4ask = decision(z_4ask, Ns(1), nb_bits, 'ask', 4, true);
    bits_decides_qpsk = decision(z_qpsk, Ns(2), nb_bits, 'psk', 4, true);
    bits_decides_8psk = decision(z_8psk, Ns(3), nb_bits, 'psk', 8, true);
    bits_decides_16qam = decision(z_16qam, Ns(4), nb_bits, 'qam', 16, true);   
    teb_bruit(:,i+1) = [TEB(bits_decides_4ask, bits, nb_bits); ...
        TEB(bits_decides_qpsk, bits, nb_bits);TEB(bits_decides_8psk, bits, nb_bits); ...
        TEB(bits_decides_16qam, bits, nb_bits)];
end

% Trace des constellations
scatterplot(symboles_4ask);
hold on
title('xe');
hold off
scatterplot(z_4ask(1:Ns:end));
hold on
title('z');
hold off

scatterplot(symboles_qpsk);
hold on
title('xe');
hold off
scatterplot(z_qpsk(1:Ns:end));
hold on
title('z');
hold off

scatterplot(symboles_8psk);
hold on
title('xe');
hold off
scatterplot(z_8psk(1:Ns:end));
hold on
title('z');
hold off

scatterplot(symboles_16qam);
hold on
title('xe');
hold off
scatterplot(z_16qam(1:Ns:end));
hold on
title('z');
hold off

% Trace du TEB pour differentes valeurs de bruit
figure(9); clf
plot(RSB, teb_bruit);
title('TEB pour differentes valeurs de bruit')
xlabel('RSB')
ylabel('TEB')
lg = legend('4-ASK', ...
    'QPSK', ...
    '8-PSK', ...
    '16-QAM', ...
    'Location','Best');
set(lg,'Interpreter','Latex');

% Comparaison TEB theorique et TEB simule
figure(10); clf
plot(RSB, teb_bruit(1,:), 'r', RSB, 3/2*qfunc(sqrt(12/15*10.^(RSB/10)))/2, 'r-.',...
    RSB, teb_bruit(2,:), 'b', RSB, 2*qfunc(sqrt(2*2*10.^(RSB/10))*sin(pi/4))/2, 'b--',...
    RSB, teb_bruit(3,:), 'g', RSB, 2*qfunc(sqrt(2*3*10.^(RSB/10))*sin(pi/8))/3, 'g--',...
    RSB, teb_bruit(4,:), 'k', RSB, 3*qfunc(sqrt(12/15*10.^(RSB/10)))/4, 'k--');
title('TEB theorique et simule des 4 chaines')
xlabel('RSB')
ylabel('TEB')
lg = legend('4-ASK simule', ...
    '4-ASK theorique', ...
    'QPSK simule', ...
    'QPSK theorique', ...
    '8-PSK simule', ...
    '8-PSK theorique', ...
    '16-QAM simule', ...
    '16-QAM theorique', ...
    'Location','Best');
set(lg,'Interpreter','Latex');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4.2.2 Comparaison des chaines de transmission %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 4.2.2.2 Comparaison des DSP des 4 chaines

% Calcul des DSP
DSP_xe_4ask = (1/length(xe_4ask))*abs(fft(xe_4ask,2^nextpow2(length(xe_4ask)))).^2;
DSP_xe_qpsk = (1/length(xe_qpsk))*abs(fft(xe_qpsk,2^nextpow2(length(xe_qpsk)))).^2;
DSP_xe_8psk = (1/length(xe_8psk))*abs(fft(xe_8psk,2^nextpow2(length(xe_8psk)))).^2;
DSP_xe_16qam = (1/length(xe_16qam))*abs(fft(xe_16qam,2^nextpow2(length(xe_16qam)))).^2;

% Trace des DSP
figure(11); clf
plot(linspace(-1,1,length(DSP_xe_4ask)), fftshift(DSP_xe_4ask), ...
    linspace(-1,1,length(DSP_xe_qpsk)), fftshift(DSP_xe_qpsk), ...
    linspace(-1,1,length(DSP_xe_8psk)), fftshift(DSP_xe_8psk), ...
    linspace(-1,1,length(DSP_xe_16qam)), fftshift(DSP_xe_16qam));
title('DSP des 4 chaines')
xlabel('frequences normalisees')
ylabel('amplitude')
lg = legend('4-ASK', ...
    'QPSK ', ...
    '8-PSK ', ...
    '16-QAM ', ...
    'Location','Best');
set(lg,'Interpreter','Latex');

