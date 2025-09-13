function [chnkrl, chnkrm, chnkrr, varargout] = get_boundary_curves(cstruct, xlims, delta, ai, xend, slope)
%
%  This subroutine returns the left, right, and middle chunkers
%  for discretized curves in the complexification dirichlet paper.
%
%  Input arguments:
%    cstruct - structure containing curve info
%       cstruct.type: "bump" - wiggly gaussian bump (default)
%                     "gotham" - Fourier series that looks like gotham
%       cstruct.disc: "unif" - use chunkerfuncuni (default)
%                     "adap" - use chunkerfunc
%       cstruct.nchs: number of chunks to use on all three segments, 
%                     must be provided if using chunkerfuncuni
%                     nchs(1) will be the number of chunks on the real
%                     part, and nchs(2) will be the number of chunks
%                     on the imaginary part
%       ctstruct.maxchunklen: maximum chunk length must be provided
%                              if using chunkerfunc
%
%   xlims: limits of real part of the curve, the real part of the curve
%         is contained in (xlims(1), xlims(2))
%   delta: parameter determining support of bump function
%          bump function is supported between (xlims(1) + delta, xlims(2) - delta)
%   ai: ramp up for imaginary part which determines its centering
%   xend: additional part in complexified parts of the curve
%         the complex parts to be discretized are (xlims(2), xlims(2) +
%         xend), and (xlims(1) - xend, xlims(1))
%   slope: slope of complexified part of curve (default (1))
%
%
    type = 'bump';
    if isfield (cstruct, 'type')
        type = cstruct.type;
    end
    disc = 'unif';
    if isfield (cstruct, 'disc')
        disc = cstruct.disc;
    end

    % set up interface
    if strcmpi(type, 'flat')
        fun1 = @(t) gaus_cos(t,0,2,1,0);
    elseif strcmpi(type, 'bump')
        fun1 = @(t) gaus_cos(t,1,4,4,0);
    elseif strcmpi(type, 'gotham')
        as = 2*[1,-0.5,0.3];
        bs = [2*pi,pi/1.3,sqrt(2)*pi];
        cs = [0.8,1.5,3.2];

        fun1 = @(t) cos_wig(t,as,bs,cs);

    end


    % set up bump function
    a0 = 2;
    b0 = 1;
    t0 = xlims(1) + delta;
    t1 = xlims(2) - delta;
    toff = sqrt(-log(eps))/a0;
    tt0 = t0 + toff;
    tt1 = t1 - toff;
    fun2 = @(t) smooth_bump(t,a0,b0,tt0,tt1);

    % set up parameters for the imaginary part
    if nargin < 5
        slope = 1;
    end
    bi = slope/6; 
    toff = 1.5*sqrt(-log(eps))/ai;

    t0i = xlims(1) - toff;
    t1i = xlims(2) + toff;
    
    fr = @(t) prod_curve(t, fun1, fun2);

    f = @(t) prod_curve(t, fun1, fun2, ai, bi, t0i, t1i);

% Start the discretization process   
    if strcmpi(disc , 'unif')
        nchs = cstruct.nchs;
        cparams.ta = xlims(1);
        cparams.tb = xlims(2);
        chnkrm = chunkerfuncuni(fr, nchs(1), cparams);
   
        cparams.ta = xlims(1) - xend;
        cparams.tb = xlims(1);
        cparams.ifclosed = 0;
        chnkrl = chunkerfuncuni(f, nchs(2), cparams);

        cparams.ta = xlims(2);
        cparams.tb = xlims(2) + xend;
        chnkrr = chunkerfuncuni(f, nchs(2), cparams);
    elseif strcmpi(disc, 'adap')
        cparams.ta = xlims(1);
        cparams.tb = xlims(2);
        cparams.eps = 1e-12;
        cparams.maxchunklen = cstruct.maxchunklen;
        chnkrm = chunkerfunc(fr,cparams);
        


        cparams.ta = xlims(1) - xend;
        cparams.tb = xlims(1);
        cparams.ifclosed = 0;
        % cparams.maxchunklen = 0.05;
        chnkrl = chunkerfunc(f,cparams);

        cparams.ta = xlims(2);
        cparams.tb = xlims(2) + xend;
        % cparams.maxchunklen = 0.05;
        chnkrr = chunkerfunc(f,cparams);
    end
    chnkrm = sort(chnkrm);
    chnkrl = sort(chnkrl);
    chnkrr = sort(chnkrr);

    varargout{1} = f;



        
end
