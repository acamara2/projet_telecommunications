function [z, xe, h] = emission_reception(symboles, Ns, alpha, span)
    % Surechantillonage 
    d_temp = kron(symboles, [1 zeros(1,Ns-1)]);

    % Introduction de zeros en debut et fin de chaine pour anticiper le retard
    % du au filtre de mise en forme
    d = zeros(1, length(d_temp)+span*Ns);

    d(span/2*Ns+1:end-span/2*Ns) = d_temp;

    % Filtrage de mise en forme 
    h = rcosdesign(alpha, span, Ns, 'sqrt');

    xe_retard = conv(d, h);

    % Elimination du retard
    xe = xe_retard(span*Ns+1:end-span*Ns);

    % Filtrage de reception
    z = filtrage_reception(xe, span, Ns, h);
end