function [grids, indexes] = setGridsAndIndices(dims)
% SETGRIDSANDINDICES Creates grids and index matrices for model computations.
%
%   INPUT:
%       dims    - Structure with model dimension parameters:
%                   .Na, .na : Coarse and fine wealth grid sizes
%                   .K       : Productivity states
%                   .B       : Amenity states
%                   .N       : Number of locations
%
%   OUTPUT:
%       grids   - Struct with state variable grids:
%                   .agrid   : Coarse asset grid (non-uniform)
%                   .ahgrid  : Fine asset grid (uniform)
%                   .eta     : Productivity grid
%                   .psi     : Relative amenity grid
%
%       indexes - Struct with index matrices and sizes for reshaping and logic:
%                   .I_s, I_a, I_N              - Indices for main grid (s,a,N)
%                   .I_sp, I_ap, I_Np, I_app    - Indices for asset choice on fine grid
%                   .I_ep, I_psip               - Indices to recover eta and psi from s
%                   .I_sm, I_am, I_Nm, I_jm     - Migration index matrices
%                   .I_hm                       - Help state index
%                   .II                         - Linear index for stayers
%                   .sz, szp, szm               - Sizes of index matrices
%
%   DEPENDENCY:
%       Requires external function `nodeunif` for asset gr9id construction.
%
%   AUTHOR: Agustin Deambrosi
%   LAST REVISED: April 2025
% =========================================================================

    %% 1. Asset Grids (coarse and fine)
    lb.a		= 0;          % Lower bound of asset holdings
    ub.a		= 30;         % Upper bound of asset holdings
    ca			= 3;          % Curvature for coarse grid

    grids.agrid	= nodeunif(dims.Na, 0, (ub.a - lb.a)^(1/ca)).^ca + lb.a;
    grids.ahgrid= nodeunif(dims.na, lb.a, ub.a);

    %% 2. Productivity and Amenity Grids
    grids.eta	= linspace(0.8, 4, dims.K);
    grids.psi	= linspace(0.7, 2.2, dims.B);

    %% 3. Indexing for V (Main state grid): (S, Na, N)
    [I_s, I_a, I_N] = ndgrid(1:dims.S, 1:dims.Na, 1:dims.N);

    %% 4. Indexing for R (Post-saving grid): (S, Na, N, na)
    [I_sp, I_ap, I_Np, I_app] = ndgrid(1:dims.S, 1:dims.Na, 1:dims.N, 1:dims.na);
    I_ep		= floor((I_s - 1) / dims.B) + 1;       % Extract eta from s
    I_psip		= mod(I_s - 1, dims.B) + 1;            % Extract psi from s

    %% 5. Migration-Related Indexing: (S, Na, N, N, H)
    [I_sr, I_ar, I_ir, I_hr] = ndgrid(1:dims.S, 1:dims.Na, 1:dims.N, 1:dims.H);
    II			= sub2ind([dims.S, dims.Na, dims.N, dims.N, dims.H], ...
                         I_sr, I_ar, I_ir, I_ir, I_hr);  % Diagonal index: stayers

    %% 6. Full Migration Indexing: (S, Na, N, N, H)
    [I_sm, I_am, I_Nm, I_jm, I_hm] = ndgrid(1:dims.S, 1:dims.Na, 1:dims.N, 1:dims.N, 1:dims.H);

    %% 7. Sizes of Indexing Structures
    sz			= size(I_s);
    szp			= size(I_sp);
    szm			= size(I_sm);

    %% 8. Pack Outputs
    indexes.I_s		= I_s;
    indexes.I_a		= I_a;
    indexes.I_N		= I_N;

    indexes.I_sp	= I_sp;
    indexes.I_ap	= I_ap;
    indexes.I_Np	= I_Np;
    indexes.I_app	= I_app;
    indexes.I_ep	= I_ep;
    indexes.I_psip	= I_psip;

    indexes.I_sm	= I_sm;
    indexes.I_am	= I_am;
    indexes.I_Nm	= I_Nm;
    indexes.I_jm	= I_jm;
    indexes.I_hm	= I_hm;
    indexes.II		= II;

    indexes.sz		= sz;
    indexes.szp		= szp;
    indexes.szm		= szm;

end
