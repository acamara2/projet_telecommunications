function bits_decides = decision(z, Ns, nb_bits, modulation, M, use_mat_func)
    % Echantillonage
    zm = z(1:Ns:end);
    
    bits_decides = zeros(1,nb_bits);
    if strcmp('psk', modulation) && M == 4 && ~use_mat_func   
        azm = real(zm);
        bzm = imag(zm);

        % Decision
        a_decides = sign(azm);
        b_decides = sign(bzm);

        % Demapping
        bits_decides(1:2:end) = a_decides;
        bits_decides(2:2:end) = b_decides;
        bits_decides = (bits_decides+1)/2;
    else
        if use_mat_func
            switch modulation
                case 'psk'
                    symboles_demod = pskdemod(zm, M, pi/4,'gray');
                case 'qam'
                    symboles_demod = qamdemod(zm,16,'gray');
                case 'ask'
                    if M == 4
                        symboles_demod = (zm >= 2)*3 + (zm >=0 & zm < 2) - (zm < -2)*3 - (zm < 0 & zm >= -2);
                        symboles_demod = symboles_demod + 2*(symboles_demod == 1) - 2*(symboles_demod == 3);
                        symboles_demod = (symboles_demod + 3)/2; 
                    else
                        error('Version M-ask uniquement disponible pour M = 4');
                    end
                otherwise
                    error('Modulation inconnue');
            end
            for i=1:log2(M)
                bits_decides(log2(M)+1-i:log2(M):end) = fix(mod(symboles_demod, 2^i)/2^(i-1));
            end
        else
            error('Version sans fonction matlab uniquement disponible pour qpsk');
        end
    end
    
    

end