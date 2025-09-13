function [f,fd,fdd] = cos_wig(t,as,bs,cs)

    f   = zeros(size(t));
    fd  = zeros(size(t));
    fdd = zeros(size(t));

    for ii=1:numel(as)
        f  = f  + as(ii)*cos(bs(ii)*t+cs(ii));
        fd = fd - as(ii)*bs(ii)*sin(bs(ii)*t+cs(ii));
        fdd= fdd- as(ii)*(bs(ii))^2*cos(bs(ii)*t+cs(ii));
    end

end