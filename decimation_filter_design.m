% Multistage CIC decimation filter design by exhaustive search.
%
% Decimating by R in one step needs a long, expensive filter. Splitting it
% into a cascade of stages, one per prime factor of R, is far cheaper --
% but the cost depends both on the order the factors are applied in and on
% how many comb sections each stage gets.
%
% This searches every ordering of the prime factors against every
% allocation of comb sections, and keeps the cheapest design whose worst
% passband deviation stays within budget.

clear all; close all; clc;

R    = 32;          % total decimation factor
yd   = 90;          % deviation budget; the constraint is delta <= 1/yd
wp   = 0.5*pi/R;    % passband edge
Kmax = 10;          % most comb sections allowed in any one stage

Rprime = factor(R);          % prime factors of R, one per stage
M      = length(Rprime);     % number of stages

% Every distinct ordering of the stages. For R = 32 the factors are all 2,
% so there is exactly one -- the ordering search only does work when R has
% distinct prime factors.
P      = unique(perms(Rprime), 'rows');
Np     = size(P, 1);

% Every allocation of comb sections across the M stages
C      = nchoosek(0:Kmax, M);
Nc     = size(C, 1);         % rows, not length(): length() returns the
                             % largest dimension, which is only the row
                             % count by coincidence here

% Aliasing bands: the frequencies that fold into the passband on
% decimation, which are the only places the response has to be small.
points_per_segment = 16;
dw   = 0.5/R/points_per_segment;
wtmp = (0:dw:1)*pi;
wq   = [];
for a = 1:length(wtmp)
    for b = 1:R/2
        if wtmp(a) > 2*b*pi/R - wp && wtmp(a) < 2*b*pi/R + wp
            wq(end+1) = wtmp(a); %#ok<SAGROW>
        end
    end
end
Q = length(wq);

deltaopt = 1/yd;             % best deviation found so far, seeded at budget
Aopt     = inf;              % operation count of the best design so far
K        = [];
Rm       = [];

for p = 1:Np
    Rp = P(p, :);

    % Operation count for each candidate allocation
    A = zeros(Nc, 1);
    for c = 1:Nc
        Kc = C(c, :);
        O  = 0;
        for m = 1:M-1
            I = prod(Rp(m+1:M));
            O = O + abs(Kc(m) - Kc(m+1)) * I;
        end
        A(c) = Kc(1)*R + O + Kc(M);
    end

    % Search cheapest-first so the break below is valid. Sort a local copy:
    % reordering C itself would leave the next ordering searching a table
    % that no longer lines up with the costs computed for it.
    [Asorted, idx] = sort(A);
    Csorted = C(idx, :);

    for c = 1:Nc
        if Asorted(c) > Aopt
            break;   % everything further down costs at least as much
        end

        delta = max(abs(freqz(Csorted(c, :), Rp, wq)));

        if delta <= deltaopt
            deltaopt = delta;
            Aopt     = Asorted(c);
            K        = Csorted(c, :);
            Rm       = Rp;
        end
    end
end

% Differential comb lengths per stage
Lk      = zeros(M, 1);
Lk(1)   = -K(1);
Lk(M)   = K(M);
for k = 2:M-1
    Lk(k) = K(k) - K(k+1);
end

fprintf('Decimation factor R = %d\n', R);
fprintf('Stage ordering:      %s\n', mat2str(Rm));
fprintf('Comb sections K:     %s\n', mat2str(K));
fprintf('Operations A:        %d\n', Aopt);
fprintf('Worst deviation:     %.6g (budget %.6g)\n', deltaopt, 1/yd);
fprintf('Aliasing points Q:   %d\n', Q);
