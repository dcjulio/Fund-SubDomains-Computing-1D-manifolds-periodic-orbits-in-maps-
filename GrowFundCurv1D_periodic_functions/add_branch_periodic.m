function manif=add_branch_periodic(manif, opts, newbranch,eigenval_idx)

if nargin < 4
    eigenval_idx = [];
end

    %--- Information of the system
    opts.thesystem       = manif.grow_info.thesystem; % What is the name of the system file
    opts.par             = manif.par; % The parameter values and names (has to match with the names defined in StdHenon3D)
    opts.user_arclength  = manif.grow_info.user_arclength; % What is the approximate arclength of the manifold
    opts.per_orbit.name  = manif.per_orbit.name;
    opts.per_orbit.coord_original = manif.per_orbit.coord_original;
    opts.stability       = manif.stability;

    %--- Accuracy parameters (default)
    %opts.accpar.alphamax=0.3;
    %opts.accpar.deltalphamax=0.001; 
    %opts.accpar.deltamin=0.000001;
    %opts.accpar.deltamax=0.01;  

    %--- Initial step (default)
    %opts.accpar.init_step=10^-7;


    % if we one one branch, we have to have at least the other branch
    if strcmp(newbranch,'neg') && isfield(manif.points{1},'pos')
        kept_branch =  'pos';
    elseif strcmp(newbranch,'pos') && isfield(manif.points{1},'neg')
        kept_branch =  'neg';
    else
        fprintf('Something went wrong... \nPlease review the only branch currently present in the manifold and ensure you are not attempting to compute it again. \n');
        return;
    end

    if strcmp(manif.orientability,'orientation-reversing')
        fprintf('The manifold is orientation-reversing and both branches are computed in the previous run\n');
        return;
    end


    fprintf('\nComputing now the %s branch...\n', newbranch);
    
    PO_idxs = 1:numel(manif.points);

    opts.branch = newbranch;
    if isempty(eigenval_idx)
        manif_new   = GrowFundCurv1D_periodic(opts);
    else
        manif_new   = GrowFundCurv1D_periodic(opts,eigenval_idx);
    end
    
    for PO_idx = PO_idxs
        manif.points{PO_idx}.(newbranch)=manif_new.points{PO_idx}.(newbranch);
    end

    % rewrite growinfo
    grow_info.(kept_branch) = manif.grow_info;     %saving old run
    grow_info.(newbranch)  = manif_new.grow_info; %saving new run
    manif.grow_info = [];
    manif.grow_info = grow_info;
    


end