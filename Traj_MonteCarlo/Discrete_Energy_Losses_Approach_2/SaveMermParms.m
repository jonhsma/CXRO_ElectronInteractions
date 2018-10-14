function SaveMermParms(fname,physparms,mermparms_lsqfit,mermparms_legend)
%%% Usage: SaveMermParms(fname,physparms,mermparms_lsqfit,mermparms_legend)
%%% fname: (INPUT) Full filename (.mat), including full path
%%% physparms: (INPUT) Physical parameters structure
%%%            ['rho','Mw','Ef','q','hbar'] [g/cm3, g/mol, eV, /m, J.s]
%%% mermparms_lsqfit: (INPUT) Oscillators in rows, each row containing 4
%%%                   columns [Ai, Ei, gamma_i]
%%% mermparms_legend: (INPUT) legend '[Ai Ei gamma_i]'

if nargin<4
    error('SaveMermParms: Need %d inputs, provided only %d',4,nargin);
end

if ~isfield(physparms,'rho')
    error('SaveMermParms: physparms structure does not contain ''rho'' value');
end

if ~isfield(physparms,'Mw')
    error('SaveMermParms: physparms structure does not contain ''Mw'' value');
end

if ~isfield(physparms,'Ef')
    error('SaveMermParms: physparms structure does not contain ''Ef'' value');
end

if ~isfield(physparms,'q')
    error('SaveMermParms: physparms structure does not contain ''q'' value');
end

if ~isfield(physparms,'hbar')
    error('SaveMermParms: physparms structure does not contain ''hbar'' value');
end

save(fname);

