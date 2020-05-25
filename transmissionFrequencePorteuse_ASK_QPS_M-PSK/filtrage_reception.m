function [z] = filtrage_reception(zb_temp, span, Ns, hr)
    % Anticipation du retard
    zb = zeros(1, length(zb_temp)+span*Ns);
    zb(span/2*Ns+1:end-span/2*Ns) = zb_temp;

    % Filtrage de reception
    z_retard = conv(zb,hr);

    % Elimination du retard
    z = z_retard(span*Ns+1:end-span*Ns);
end