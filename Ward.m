function [Y_ward, I_ward, Vb, Vi] = Ward(system, External, Bound, Internal, Ie, Ib, Ii)
define_constants;

bus = system.bus;
[~, Y] = My_makeY(system);

External = External(:).';
Bound = Bound(:).';
Internal = Internal(:).';

all_nodes = [External, Bound, Internal];
[~, loc] = ismember(all_nodes, bus(:, BUS_I));
Y_EBI = Y(loc, loc);

nE = numel(External);
nB = numel(Bound);
nI = numel(Internal);

iE = 1:nE;
iB = nE + (1:nB);
iI = nE + nB + (1:nI);

Y_EE = Y_EBI(iE, iE);
Y_EB = Y_EBI(iE, iB);
Y_BE = Y_EBI(iB, iE);
Y_BB = Y_EBI(iB, iB);
Y_BI = Y_EBI(iB, iI);
Y_IB = Y_EBI(iI, iB);
Y_II = Y_EBI(iI, iI);

Ie = Ie(:);
Ib = Ib(:);
Ii = Ii(:);

YEEinv_YEB = Y_EE \ Y_EB;
YEEinv_IE  = Y_EE \ Ie;

Y_BB_t = Y_BB - Y_BE * YEEinv_YEB;
I_B_t  = Ib    - Y_BE * YEEinv_IE;

Y_ward = [Y_BB_t, Y_BI; Y_IB, Y_II];
I_ward = [I_B_t; Ii];

V = Y_ward \ I_ward;
Vb = V(1:nB);
Vi = V(nB+1:end);
end
