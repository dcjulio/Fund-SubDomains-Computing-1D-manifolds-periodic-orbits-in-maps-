function Manif=GrowFundCurv1D_periodic(opts,eigenval_idx)
% Input:
% - opts: all the options needed to compute the manifold
% - eigenval_idx: optional argument. It is which eigenvector to consider for 2D manifolds
if nargin < 2
    eigenval_idx = [];
else
    if eigenval_idx~=1 && eigenval_idx~=2
        fprintf('\n ERROR! The optional second argument eigenval_idx can be either 1 or 2\n')
        return;
    end
end
    
    %% Initializing the manifold structure
    manif = init_manif_periodic(opts);
    % figure
    % hold on

    tic
    %%
    %saves the information of the system stored in the structure 'manif' and 'opts'
    par=manif.par;
    thesystem=opts.thesystem;
    per_orbit=manif.per_orbit.coord_compactified;
    period = numel(per_orbit.x);
    name = opts.per_orbit.name;

    %% Warnings in case the specifications are not correct for the computation of this manifold

    if strcmp(manif.orientability,'orientation-reversing')
        fprintf('\n----  Manifold is non-orientable. We will proceed computing both branches\n');
        mapiter=2;
    else
        mapiter=1;
    end

    % if computing unstable manifold, use image
    % if computing stable manifold, use preimage
    if strcmp(manif.stability,'Smanifold')
        sign = -1; %inverse
    elseif strcmp(manif.stability,'Umanifold')
        sign = 1; 
    end
    
    if manif.grow_info.init_step < 0
        branch1='neg';
        branch2='pos';
    else
       branch1='pos';
       branch2='neg';
    end

    %% Warnings if the manifold we want to compute is not 1D
    
    if manif.grow_info.dimension>1
        fprintf('\n----  The dimension of the manifold is %i and this algorithm is designed to compute one-dimensional manifolds. -----\n',manif.grow_info.dimension)
        fprintf(' To choose automatically which of the two eigenvalues to use, include a second argument to GrowFundCurv1D(opts,#) where #: 1 or 2\n\n') 
        fprintf('The eigenvectors are:\n')

        for i=1:manif.grow_info.dimension
            fprintf('(%i) Eigenvalue: %.3f, eigenvector of the first element of the per. orbit: (%.3f, %.3f, %.3f)\n',i, manif.grow_info.eigval(i), manif.grow_info.eigvecs(1,i,1),manif.grow_info.eigvecs(2,i,1),manif.grow_info.eigvecs(3,i,1))
        end

        if ~isempty(eigenval_idx)
            x = eigenval_idx;
            fprintf('\nChoice: eigenval_idx = %i\n\n', eigenval_idx)
        else
            fprintf('\n Which eigenvector do you want to use?') 
            prompt = "\n... Press 0 in case you want to exit \n\n";
            x = input(prompt);
        end

        if x == 0
            return
        else
            manif.grow_info.eigval=manif.grow_info.eigval(x);
            manif.grow_info.eigvecs=manif.grow_info.eigvecs(:,x,:);
        end

    end
%% Error if is a repeller or attractor
if numel(manif.grow_info.eigvecs(:,:,1))==0
    if strcmp(manif.stability,'Smanifold')
        fprintf('Error: This periodic orbit appears to be a repeller and it doesnt have a Stable manifold.\n');
    else
        fprintf('Error: This periodic orbit appears to be a attractor and it doesnt have an Unstable manifold.\n');
    end
    return
end

%% Printing the general information of this run

    names=fieldnames(par);
    
    fprintf('\n----');
    for k=1:length(names)
        fprintf(' %s:%0.2f ',names{k},par.(names{k}));
    end
    fprintf('----');

    if strcmp(manif.orientability,'orientation-preserving')
        fprintf('\n----  Computing W%s(%s), %s branch up to arclength %.0f each ---- ',lower(opts.stability(1)),opts.per_orbit.name, branch1, opts.user_arclength);
    else
        fprintf('\n----  Computing W%s(%s), both branches up to arclength %.0f each ---- ',lower(opts.stability(1)),opts.per_orbit.name,opts.user_arclength);
    end
    
    fprintf('\n----  Total maximum iterations of the fundamental domain: %i ---- ',manif.grow_info.max_iter);
    fprintf('\n\n')
    fprintf('----- Acc. Conditions ----- ')
    fprintf('\n| AlphaMax:%.2e',manif.grow_info.alphamax)
    fprintf('\n| DeltaAlphaMax:%.2e ',manif.grow_info.deltalphamax)
    fprintf('\n| DeltaMin:%.2e ' ,manif.grow_info.deltamin)
    fprintf('\n| DeltaMax:%.2e \n' ,manif.grow_info.deltamax)
    fprintf('---------- \n')
    
    
    %-------------------------------------------------------------------
    %-------------------------------------------------------------------
    % initializing the fields for the initial information
    manif.grow_info.runinf.rem_deltamin=0; %how many points are removed because of deltamin
    manif.grow_info.runinf.rem_nan=0; %how many points are removed because of NaN values
    manif.grow_info.runinf.rem_inf=0; %how many points are removed at infinity because of duplication
    manif.grow_info.runinf.add_alphamax=0; %how many points are added because of alpha_max
    manif.grow_info.runinf.add_deltamax=0;%how many points are added because of delta_max
    manif.grow_info.runinf.add_deltalphamax=0;  %how many points are added because of (delta alpha)_max
%% Mapping rules

    % Iteration order of periodic orbits (naturally for Unstable manifolds)
    if strcmp(manif.stability,'Umanifold')
        order_periodic    = 1:period; %[1 2 3... k]
    elseif strcmp(manif.stability,'Smanifold') % if we are computing Stable, we use the inverse
        order_periodic = [1 k:-1:2]; %[1 k k-1 ... 2]
    end
    preimage_periodic = circshift(order_periodic, 1); %shift right (last becomes first)
    % [k, 1, 2, 3,   ..., k-1] Unstable manif
    % [2, 1, k, k-1, ..., 3] Stable manif

%% computes first fundamental domain


% by default, we choose the eigenvector going to positive x for the first
% element of the per orbit
if manif.grow_info.eigvecs(1,:,1) < 0 % if it is going to negative x, then consider the other side
    manif.grow_info.eigvecs = -manif.grow_info.eigvecs;
end

% Choosing the initial distance correctly:
% To avoid problems when the contraction/expansion rate is too strong/weak,
% We define an starting segments from each periodic point up to
% opts.accpar.init_step in the direction of each eigenvector
% Then, we iterate those starting pieces to define the first fundamental domain.

init_step = manif.grow_info.init_step;

% Initialize structure of the stored fundamental domain
fdom.points = cell(1,period);
for K = 1:period
    k = order_periodic(K);
    step_x = init_step*manif.grow_info.eigvecs(1,:,k); %(x,eigv,per)
    step_y = init_step*manif.grow_info.eigvecs(2,:,k); %(y,eigv,per)
    step_z = init_step*manif.grow_info.eigvecs(3,:,k); %(z,eigv,per)

    % interpolate N=4 points linearly
    fdom.points{k}.(branch1).x = linspace(per_orbit.x(k), per_orbit.x(k)+step_x, 6);
    fdom.points{k}.(branch1).y = linspace(per_orbit.y(k), per_orbit.y(k)+step_y, 6);
    fdom.points{k}.(branch1).z = linspace(per_orbit.z(k), per_orbit.z(k)+step_z, 6);

    fdom.points{k}.(branch1).arc = arclength(fdom.points{k}.(branch1));
    fdom.points{k}.(branch1).idx_fdom = [1 numel(fdom.points{k}.(branch1).x)];

    if strcmp(manif.orientability,'orientation-reversing')
        % interpolate N=4 points linearly
        fdom.points{k}.(branch2).x = linspace(per_orbit.x(k), per_orbit.x(k)-step_x, 6);
        fdom.points{k}.(branch2).y = linspace(per_orbit.y(k), per_orbit.y(k)-step_y, 6);
        fdom.points{k}.(branch2).z = linspace(per_orbit.z(k), per_orbit.z(k)-step_z, 6);
    
        fdom.points{k}.(branch2).arc = arclength(fdom.points{k}.(branch2));
        fdom.points{k}.(branch2).idx_fdom = [1 numel(fdom.points{k}.(branch2).x)];
    end

end

%% ///////////////// FIG
% fig_handles = gobjects(1, period);
% 
% for K = 1:period
%     k = order_periodic(K);
%     name_plot = manif.per_orbit.name+"_"+k;
%     fig_handles(k) = figure('Name', name_plot);
%     hold on;
% 
%     set(groot, 'CurrentFigure', fig_handles(k));
%     plot3(per_orbit.x(k),per_orbit.y(k),per_orbit.z(k),'*');
%     hold on
%     plot3(fdom.points{k}.(branch1).x,fdom.points{k}.(branch1).y,fdom.points{k}.(branch1).z,'.-','MarkerSize',15, 'LineWidth',1.5);
%     if strcmp(manif.orientability,'orientation-reversing')
%         plot3(fdom.points{k}.(branch2).x,fdom.points{k}.(branch2).y,fdom.points{k}.(branch2).z,'.-','MarkerSize',15, 'LineWidth',1.5);
%     end
%     text(manif.per_orbit.coord_compactified.x(k),manif.per_orbit.coord_compactified.y(k),manif.per_orbit.coord_compactified.z(k),manif.per_orbit.name+"\_"+k);
%     xlabel('x')
%     ylabel('y')
%     zlabel('z')
% end
%/////////////////
%%


    % Initializing the structures
    manif.points      = cell(1,period);
    total_arc_branch1 = zeros(1,period);
    total_arc_branch2 = zeros(1,period);

    [manif.points{1:period}] = deal(struct('name',[], branch1, [])); % initialize the cells
    for K = 1:period
        k = order_periodic(K);
        manif.points{k}.name = sprintf('%s_%s',name,num2str(k));

        manif.points{k}.(branch1).x   = fdom.points{k}.(branch1).x;
        manif.points{k}.(branch1).y   = fdom.points{k}.(branch1).y;
        manif.points{k}.(branch1).z   = fdom.points{k}.(branch1).z;
        manif.points{k}.(branch1).arc = fdom.points{k}.(branch1).arc;
        manif.points{k}.(branch1).idx_fdom = fdom.points{k}.(branch1).idx_fdom;

        manif.points{k}.(branch1).branch_preimage = branch1; %name of the branch
        manif.points{k}.(branch1).name_preimage   = [];
        manif.points{k}.(branch1).idx_preimages   = [];

        total_arc_branch1(k) =  manif.points{k}.(branch1).arc(end);
    end
    for K = 1:period
        k = order_periodic(K);
        k_pre = preimage_periodic(K);
        manif.points{k}.(branch1).name_preimage = manif.points{k_pre}.name;
        manif.points{k}.(branch1).idx_preimages = zeros(size(manif.points{k_pre}.(branch1).x));
    end



 if strcmp(manif.orientability,'orientation-reversing') 

    for K = 1:period
        k = order_periodic(K);

        manif.points{k}.(branch2).x=fdom.points{k}.(branch2).x;
        manif.points{k}.(branch2).y=fdom.points{k}.(branch2).y;
        manif.points{k}.(branch2).z=fdom.points{k}.(branch2).z;
        manif.points{k}.(branch2).arc=fdom.points{k}.(branch2).arc;
        manif.points{k}.(branch2).idx_fdom=fdom.points{k}.(branch2).idx_fdom;

        manif.points{k}.(branch2).branch_preimage = []; %name of the branch
        manif.points{k}.(branch2).name_preimage   = [];
        manif.points{k}.(branch2).idx_preimages   = [];

        total_arc_branch2(k) =  manif.points{k}.(branch2).arc(end);

    end
    for K = 1:period
        k = order_periodic(K);
        k_pre = preimage_periodic(K);

        manif.points{k}.(branch2).name_preimage = manif.points{k_pre}.name;
        if k == 1
            manif.points{k}.(branch1).branch_preimage = branch2;
            manif.points{k}.(branch2).branch_preimage = branch1;
        else
            manif.points{k}.(branch1).branch_preimage = branch1;
            manif.points{k}.(branch2).branch_preimage = branch2;
        end
    end

 end

Manif = manif;
% 
% 
% progress table
varNames = arrayfun(@(i) sprintf('p%d_%d', period, i), order_periodic, 'UniformOutput', false);
if strcmp(manif.orientability,'orientation-preserving')
    data = num2cell(zeros(1, period), 1);
    T = table(data{:},'VariableNames',varNames,'RowName',{sprintf('%s arc%f\n',branch1)});
    T{1,:} = round(total_arc_branch1, 2, 'significant');
else
    data = num2cell(zeros(2, period), 1);
    T = table(data{:},'VariableNames',varNames,'RowName',{sprintf('%s arc%f\n',branch1),sprintf('%s arc%f\n',branch2)});
    T{1,:} = round(total_arc_branch1, 2, 'significant');
    T{2,:} = round(total_arc_branch2, 2, 'significant');
end
disp(T)
% 
% 
% 
%----------- Starting the loop
% warning = 1; %to know if we went to some warnning
iter = 0;
stop_convergence_branch1 = zeros(1,period); % flag to stop iterating if we branch1 reached convergence to something
stop_convergence_branch2 = zeros(1,period); % flag to stop iterating if we branch2 reached convergence to something
stop_arc = 0; % flag manifold has reached the arclength
warning=0;

% iter: number of iterations
while (iter < manif.grow_info.max_iter) && (stop_arc == 0) %&& ((sum(stop_convergence_branch1) + sum(stop_convergence_branch2)) ~= mapiter*period)
       % we haven't reached the max number of iterations
       % and, we haven't reached the arclength (orientation-preserving case)
       % and, we haven't reached the arclength for the two branches of the manifolds (orientation-reversing case)

    iter = iter + 1;
    periodic_idx      = order_periodic(mod(iter, period)+1); %index if current segment of manifold. It is which periodic point we are computing now
    preimage_periodic_idx = preimage_periodic(mod(iter, period)+1); % where this periodic point comes from

    % for stable manifolds we go backwards, for unstable, forward

    %Chech in which branch are we at (for orientation reversing, we first
    %compute branch 1 and then branch 2)
    if strcmp(manif.orientability,'orientation-preserving')
        branch = branch1;
    elseif strcmp(manif.orientability,'orientation-reversing')
        if mod(floor(iter/period), 2) == 0 % 0: branch1, 1:branch2
            branch = branch1;
        else
            branch = branch2;
        end
    end

    % we take the image of the last fdom of the preimage branch
    fdom_level = floor((iter - 1)/(mapiter*period)) + 1;

    branch_preimage = Manif.points{periodic_idx}.(branch).branch_preimage;
    pre_fdom1 = Manif.points{preimage_periodic_idx}.(branch_preimage).idx_fdom(fdom_level,1);
    pre_fdom2 = Manif.points{preimage_periodic_idx}.(branch_preimage).idx_fdom(fdom_level,2);

    fdom_points.x = Manif.points{preimage_periodic_idx}.(branch_preimage).x(pre_fdom1:pre_fdom2);
    fdom_points.y = Manif.points{preimage_periodic_idx}.(branch_preimage).y(pre_fdom1:pre_fdom2);
    fdom_points.z = Manif.points{preimage_periodic_idx}.(branch_preimage).z(pre_fdom1:pre_fdom2);


    % mapping the points
    mappoints = thesystem.mapping(fdom_points,opts,sign);
    idx_preimages   = pre_fdom1:pre_fdom2;

%/////////////////
% if iter<=mapiter*period
%     set(groot, 'CurrentFigure', fig_handles(periodic_idx));
%     plot3(mappoints.x,mappoints.y,mappoints.z,'ko-', 'MarkerSize',10);
% else
%     set(groot, 'CurrentFigure', fig_handles(periodic_idx));
%     plot3(mappoints.x,mappoints.y,mappoints.z,'k*-', 'MarkerSize',10);
% end

%/////////////////


     %% STARTING THE ALGORITHM

    %% ----------- Removing points at infinity or NaN values. This is only because issues from compactification

    %Delete NaN values
    nan_idx=union(union(find(isnan(mappoints.x)),find(isnan(mappoints.y))),find(isnan(mappoints.z)));
    Manif.grow_info.runinf.rem_nan=Manif.grow_info.runinf.rem_nan+numel(nan_idx);

      if numel(nan_idx)==numel(mappoints.x)
        fprintf('\n\n ALL THE POINTS HAVE BEEN MAPPED TO NAN VALUES. LAST POINT IS AT (%f,%f,%f)\n\n',fdom_points.x(end),fdom_points.y(end),fdom_points.z(end))

        break;
       end

    mappoints.x(nan_idx)=[]; 
    mappoints.y(nan_idx)=[]; 
    mappoints.z(nan_idx)=[]; 

    fdom_points.x(nan_idx)=[]; 
    fdom_points.y(nan_idx)=[]; 
    fdom_points.z(nan_idx)=[]; 

    idx_preimages(nan_idx)=[]; 


    %Delete points at infinity
    inf_idx=union(find((mappoints.x.^2+mappoints.y.^2)==1),find(abs(mappoints.z)==1));
    Manif.grow_info.runinf.rem_inf=Manif.grow_info.runinf.rem_inf+numel(inf_idx);


    if numel(inf_idx) == numel(mappoints.x)
        fdom_points.x = fdom_points.x(end);
        fdom_points.y = fdom_points.y(end);
        fdom_points.z = fdom_points.y(end);

        mappoints.x = mappoints.x(end);
        mappoints.y = mappoints.y(end);
        mappoints.z = mappoints.y(end);

        idx_preimages=idx_preimages(1);
    else

        mappoints.x(inf_idx)=[]; 
        mappoints.y(inf_idx)=[]; 
        mappoints.z(inf_idx)=[]; 
    
        fdom_points.x(inf_idx)=[]; 
        fdom_points.y(inf_idx)=[]; 
        fdom_points.z(inf_idx)=[]; 
    
        idx_preimages(inf_idx)=[];
    end


     %% Replace first point of the current mapped segment by the last point of the previous segment (continuous manifold)
     % and monitoring the distance between first point of the mapped points and the last point of previous fundamental domain on that branch



    if iter > period*mapiter % it is not the first time to store an fdom after the initial segment 

        % concatenate the new fdom with the previous fdom

        %first check the distance between the first point of the new fdom
        %with the last point of the previous fdom
        fdom_step = sqrt((Manif.points{periodic_idx}.(branch).x(end) - mappoints.x(1))^2 +...
        (Manif.points{periodic_idx}.(branch).y(end) - mappoints.y(1))^2 +...
        (Manif.points{periodic_idx}.(branch).z(end) - mappoints.z(1))^2);
        % fprintf('dist %.e\n', fdom_step)


        % check distance
        if fdom_step > Manif.grow_info.deltamin && warning < 3
            fprintf('Warning! The distance between the first point of the current fdom and the last point of the previous fdom exceeds Deltamin. Current distance: %.e %.e\n', fdom_step)
            prompt = "\n... Press sany key to continue, or 0 to stop receiving this message\n\n";
            x = input(prompt);
            warning = warning +1;
        end

         %replace the first point of the mapped points by the last point of the previous segment
         mappoints.x(1) = Manif.points{periodic_idx}.(branch).x(end);
         mappoints.y(1) = Manif.points{periodic_idx}.(branch).y(end);
         mappoints.z(1) = Manif.points{periodic_idx}.(branch).z(end);

     end


    %%
    fprintf('\n FUNDAMENTAL SUBDOMAIN %i of %s %s branch',fdom_level,manif.points{periodic_idx}.name, branch)
    fprintf(' (%i total iterations) \n', iter);


    arc_mappoints = arclength(mappoints); %arclength of the mapped fundamental domain

    % Check how much arclength is left to compute
    needed_arc = opts.user_arclength - (sum(total_arc_branch1) + sum(total_arc_branch2));

% 
%     % Check how long is the new fundamental domain, if is less than
%     % Deltamin, then mark it as convergence.
%     if arc_mappoints(end) < abs(manif.grow_info.init_step) && iter>3*mapiter*period
% % 
%         fprintf(' NOTE: The arclength of this fundamental domain is less than initial step |delta|... \n');
% stop_convergence_branch1
%         if strcmp(manif.orientability,'orientation-preserving')
%             stop_convergence_branch1(periodic_idx) = 1;
%         elseif strcmp(manif.orientability,'orientation-reversing')
%             if mod(floor(iter/period), 2) == 0 % 0: branch1, 1:branch2
%                 stop_convergence_branch1(idx_seg) = 1;
%             else
%                 stop_convergence_branch2(idx_seg)=1;
%             end
    %     end
    % end
% 

    % if in this iteration we exceed the desired arclength, then we chop
    % the fdomain up to the desired arclength
    if arc_mappoints(end) > needed_arc

         idx_arc = find(arc_mappoints > needed_arc,1); %where we exceed the extra needed arc
         % chop the fund domain and the mappoints up to there

         fdom_points.x = fdom_points.x(1:idx_arc);
         fdom_points.y = fdom_points.y(1:idx_arc);
         fdom_points.z = fdom_points.z(1:idx_arc);

         mappoints.x = mappoints.x(1:idx_arc);
         mappoints.y = mappoints.y(1:idx_arc);
         mappoints.z = mappoints.z(1:idx_arc);

         idx_preimages = idx_preimages(1:idx_arc);

         stop_arc = 1; % The manifold has reached the arclength and this is the last iteration
    end
% 
% 
%------------------------------------------------------- 
%---%----------- Adding points depending on Acc. Cond.  

    %initializing the structures to add points
    add_acc = struct(); 
    newpoints = struct();
    mapnewpoints = struct();


	add_acc.iter   = 0;
    add_acc.failed = []; % points that failed acc cond in last loop
    add_acc.loop   = true; % still doing the while loop 

    %if is the NOT the first fdom on that branch
    %we also take into account the point at the end of the previous fdom 
    if iter > period*mapiter
        mappoints.x = [Manif.points{periodic_idx}.(branch).x(end-1) mappoints.x];
        mappoints.y = [Manif.points{periodic_idx}.(branch).y(end-1) mappoints.y]; 
        mappoints.z = [Manif.points{periodic_idx}.(branch).z(end-1) mappoints.z]; 

        % we add a dummy at the beggining, we NEVER interpolate between fdom_points.x(1) and fdom_points.x
        fdom_points.x = [fdom_points.x(1) fdom_points.x];
        fdom_points.y = [fdom_points.y(1) fdom_points.y]; 
        fdom_points.z = [fdom_points.z(1) fdom_points.z]; 
    end


    t_initial = 0:1/(numel(fdom_points.x)-1):1; % parametrization for meshpoints

%---%-------------- Loop of the same mesh checking acc cond (this adds points)
    while add_acc.loop 
%---%--------------
        add_acc.loop = false; %to stop while loop % if at least one point is added this turns true
        add_acc.iter = add_acc.iter+1; %counter of iterations


        % Interpolating points from previous fundamental domain.
         if add_acc.iter == 1
            tt        = t_initial;
            t_interp  = tt(1:end-1)+(tt(2:end)-tt(1:end-1))/2; % parametrization of interpolated points
            interp    = makima3D(fdom_points,t_initial,t_interp); % compute interpolated preimage
            mapinterp = thesystem.mapping(interp,opts,sign); % interpolated image
         else
             tt        = sort([tt t_interp(add_acc.add)]); %parametrization of (new) mesh points
             t_interp  = tt(1:end-1)+(tt(2:end)-tt(1:end-1))/2; % parametrization of (new) interpolated points
             interp    = makima3D(fdom_points,t_initial,t_interp); % compute interpolated preimage
             mapinterp = thesystem.mapping(interp,opts,sign); % interpolated image
         end

%%
        add_acc.add = []; %points we are going to add

        % idx of points to check acc cond
        if add_acc.iter == 1
            for_idx = 2:(numel(mappoints.x)-1);
        else
            for_idx = add_acc.failed;
        end      

        fprintf('  loop number %i (points to check %i...',add_acc.iter,numel(for_idx));

%-------%---------- Going through the points that failed
        million=0;
        for k = for_idx
%-------%----------      
            % a flag for when # million points have been checked
            million=million+1;

            if floor(million/1000000)==ceil(million/1000000) %is integer?
                fprintf(' -checkpoint %i million points checked-...',floor(million/1000000));
            end


            % coordinates of mapped points
            add_acc.p0 = [mappoints.x(k-1), mappoints.y(k-1), mappoints.z(k-1)];
            add_acc.p1 = [mappoints.x(k),   mappoints.y(k),   mappoints.z(k)]; % the point we are actually looking at
            add_acc.p2 = [mappoints.x(k+1), mappoints.y(k+1), mappoints.z(k+1)];

            % Distance btw points
            add_acc.delta0 = norm(add_acc.p1-add_acc.p0); % before
            add_acc.delta2 = norm(add_acc.p1-add_acc.p2); % after
            add_acc.alpha  = angles(add_acc.p0,add_acc.p1,add_acc.p2); % angle btw points

            % points btw p0p1 and p1p2 in the interpolated points
            add_acc.p0_new = [mapinterp.x(k-1), mapinterp.y(k-1), mapinterp.z(k-1)];
            add_acc.p2_new = [mapinterp.x(k),   mapinterp.y(k),   mapinterp.z(k)];


%-----------%------ Adding points


%-----------%------ If delta > deltamax 

            %------ If it is the second point in the mesh and first segment in the branch
            %------ Check the first delta and add a point before if needed
            if k==2 && iter < mapiter*period && add_acc.delta0>manif.grow_info.deltamax  

                %add point p01
                add_acc.add  = [add_acc.add k-1]; %idx of the point we are going to add
                add_acc.loop = true; % We have to check if we need to put more points in the mesh
                Manif.grow_info.runinf.add_deltamax = Manif.grow_info.runinf.add_deltamax+1;

            %------ Check the second delta and add a point after 
            elseif add_acc.delta2>manif.grow_info.deltamax %k>2 && add_acc.delta2>manif.grow_info.deltamax
                %add point p12
                add_acc.add  = [add_acc.add k]; %idx of the point we are going to add
                add_acc.loop = true; %idx of the point we are going to add
                Manif.grow_info.runinf.add_deltamax = Manif.grow_info.runinf.add_deltamax+1;



%-----------%------ If alpha > alphamax  or   Delta*alpha > Delta*alpha max

            %------ If alpha fails or BOTH Delta*alpha fail
            %------ Choose where to add a point.
            %-- Only if: we are either on the first segment, or k is not 2 (iter < mapiter*period || k~=2)
            elseif (iter < mapiter*period || k~=2) && (add_acc.alpha>=manif.grow_info.alphamax || (add_acc.delta0*add_acc.alpha>=manif.grow_info.deltalphamax && add_acc.delta2*add_acc.alpha>=manif.grow_info.deltalphamax)) %1 

            %------ Add point where Delta>Deltamin

                %-- If only Delta0>Deltamin
                %-- Add a point btw p0 and p1
                if add_acc.delta0>manif.grow_info.deltamin && add_acc.delta2<manif.grow_info.deltamin 
                    % Add a point if we didnt added the point in the previous acc cond checks
                    if numel(add_acc.add)==0 || add_acc.add(end)~=k-1
                        add_acc.add  = [add_acc.add k-1];  %idx of the point we are going to add
                        add_acc.loop = true; %idx of the point we are going to add
                        if (add_acc.delta0*add_acc.alpha>=manif.grow_info.deltalphamax && add_acc.delta2*add_acc.alpha>=manif.grow_info.deltalphamax) 
                            Manif.grow_info.runinf.add_deltalphamax = Manif.grow_info.runinf.add_deltalphamax+1;
                        else
                            Manif.grow_info.runinf.add_alphamax = Manif.grow_info.runinf.add_alphamax+1;
                        end
                    end
                end

                %-- If only Delta2>Deltamin
                %-- Add a point btw p1 and p2
                if add_acc.delta2>manif.grow_info.deltamin && add_acc.delta0<manif.grow_info.deltamin
                    %add point p12
                    add_acc.add  = [add_acc.add k]; %idx of the point we are going to add
                    add_acc.loop = true; %idx of the point we are going to add
                    if (add_acc.delta0*add_acc.alpha>=manif.grow_info.deltalphamax && add_acc.delta2*add_acc.alpha>=manif.grow_info.deltalphamax) 
                        Manif.grow_info.runinf.add_deltalphamax = Manif.grow_info.runinf.add_deltalphamax+1;
                    else
                        Manif.grow_info.runinf.add_alphamax = Manif.grow_info.runinf.add_alphamax+1;
                    end
                end

                %-- If both Delta0 and Delta 2 > Deltamin
                %-- Choose where to add point
                if add_acc.delta2>manif.grow_info.deltamin && add_acc.delta0>manif.grow_info.deltamin
                    add_acc.alpha0_new = angles(add_acc.p0_new,add_acc.p1,add_acc.p2); % angle btw points
                    add_acc.alpha2_new = angles(add_acc.p0,add_acc.p1,add_acc.p2_new); 

                    if add_acc.alpha0_new < add_acc.alpha2_new 
                        % Add a point if we didnt added the point in the previous acc cond checks
                        if numel(add_acc.add)==0 || add_acc.add(end)~=k-1
                            add_acc.add  = [add_acc.add k-1];  %idx of the point we are going to add
                            add_acc.loop = true; %idx of the point we are going to add
                            if (add_acc.delta0*add_acc.alpha>=manif.grow_info.deltalphamax && add_acc.delta2*add_acc.alpha>=manif.grow_info.deltalphamax) 
                                Manif.grow_info.runinf.add_deltalphamax = Manif.grow_info.runinf.add_deltalphamax+1;
                            else
                                Manif.grow_info.runinf.add_alphamax = Manif.grow_info.runinf.add_alphamax+1;
                            end
                        end
                    else
                        %add point p12
                        add_acc.add  = [add_acc.add k]; %idx of the point we are going to add
                        add_acc.loop = true; %idx of the point we are going to add
                        if (add_acc.delta0*add_acc.alpha>=manif.grow_info.deltalphamax && add_acc.delta2*add_acc.alpha>=manif.grow_info.deltalphamax) 
                            Manif.grow_info.runinf.add_deltalphamax = Manif.grow_info.runinf.add_deltalphamax+1;
                        else
                            Manif.grow_info.runinf.add_alphamax = Manif.grow_info.runinf.add_alphamax+1;
                        end
                    end
                end



%-----------%------ If just one Delta*alpha > Delta*alpha max ( and alpha < alphamax (previous elseif is when alpha > alphamax )

            %------ Delta0*alpha > max, and Delta0 > Deltamin
            %------ Add a point btw p0 and p1
            %-- Only if: we are either on the first segment, or k is not 2 (iter < mapiter || k~=2)
            elseif (iter < mapiter*period || k~=2) && (add_acc.delta0*add_acc.alpha>=manif.grow_info.deltalphamax && add_acc.delta0>manif.grow_info.deltamin)
                % Add a point if we didnt added the point in the previous acc cond checks
                if numel(add_acc.add)==0 || add_acc.add(end)~=k-1
                    add_acc.add  = [add_acc.add k-1];  %idx of the point we are going to add
                    add_acc.loop = true; %idx of the point we are going to add
                    Manif.grow_info.runinf.add_deltalphamax = Manif.grow_info.runinf.add_deltalphamax+1;
                end

            %------ Delta2*alpha > max, and Delta2 > Deltamin
            %------ Add a point btw p0 and p1
            elseif add_acc.delta2*add_acc.alpha>=manif.grow_info.deltalphamax && add_acc.delta2>manif.grow_info.deltamin
                %add point p12
                add_acc.add  = [add_acc.add k]; %idx of the point we are going to add
                add_acc.loop = true; %idx of the point we are going to add
                Manif.grow_info.runinf.add_deltalphamax = Manif.grow_info.runinf.add_deltalphamax+1;


%-----------%------              
            end    % (if loop) Adding points     
%-----------%------      
%-------%---------- 
        end       % (for k = for_idx) Going through the points that failed 
%-------%----------

        newpoints.x = interp.x(add_acc.add);
        newpoints.y = interp.y(add_acc.add);
        newpoints.z = interp.z(add_acc.add);

        mapnewpoints.x = mapinterp.x(add_acc.add);
        mapnewpoints.y = mapinterp.y(add_acc.add);
        mapnewpoints.z = mapinterp.z(add_acc.add);

        new_idx_preimages = idx_preimages(add_acc.add); % from where we interpolated from

       fprintf(' added points: %i) \n',numel(add_acc.add));


        % Add the points 
        if ~isempty(add_acc.add)

            % get updated idx of failed points
            plus           = 0:(numel(add_acc.add)-1);
            add_acc.failed = unique(sort([add_acc.add+plus add_acc.add+plus+1 add_acc.add+plus+2]));
            add_acc.failed = add_acc.failed(add_acc.failed>1);
            add_acc.failed = add_acc.failed(add_acc.failed<numel(mappoints.x)+numel(add_acc.add));

            % add points in the mapped manifold and in the old manifold
            mappoints.x = insert(mappoints.x,mapnewpoints.x,add_acc.add);
            mappoints.y = insert(mappoints.y,mapnewpoints.y,add_acc.add);
            mappoints.z = insert(mappoints.z,mapnewpoints.z,add_acc.add); 

            idx_preimages = insert(idx_preimages,new_idx_preimages,add_acc.add);

        end

if add_acc.iter > opts.max_refines
    add_acc.loop = false;
end

%---%--------------       
    end           % (while add_acc.loop)
%---%--------------


    % if is the NOT the first fdom on that branch, then update
    % mappoint so it doesn't contain previous fdom
    if iter > mapiter*period
        mappoints.x = mappoints.x(2:end);
        mappoints.y = mappoints.y(2:end);
        mappoints.z = mappoints.z(2:end);
    end
% 
%         %/////////////////
%         % set(groot, 'CurrentFigure', fig_handles(idx_seg));
%         % plot3(mappoints.x,mappoints.y,mappoints.z,'k.-');
%         %/////////////////
% 
    % if we only have the initial segment on the branch
    % then chop the overlaping part
    if iter <= period*mapiter
%         % %/////////////////
%         % set(groot, 'CurrentFigure', fig_handles(idx_seg));
%         % plot3(mappoints.x,mappoints.y,mappoints.z,'r.-');
%         % %/////////////////
% 
% 
%         % fprintf('\n\nchop chop chop\n\n')
%         % dist_init_segment = sqrt((Manif.points{idx_seg}.(branch).x(end) - per_orbit.x(idx_seg))^2 +...
%         %                          (Manif.points{idx_seg}.(branch).y(end) - per_orbit.y(idx_seg))^2 +...
%         %                          (Manif.points{idx_seg}.(branch).z(end) - per_orbit.z(idx_seg))^2);
        dist_init_segment = abs(Manif.grow_info.init_step);
%         % fprintf('\n\ndist init segment %e\n\n',dist_init_segment)
        dist_new_segment = sqrt((mappoints.x - per_orbit.x(periodic_idx)).^2 +...
                                (mappoints.y - per_orbit.y(periodic_idx)).^2 +...
                                (mappoints.z - per_orbit.z(periodic_idx)).^2);
% 
        idx = find(dist_new_segment > dist_init_segment,1,'first'); % find first time the dist of the new segment is greater than the last point of previous segment

        % distance from last point of previous segment to the line
        P0 = [Manif.points{periodic_idx}.(branch).x(end), Manif.points{periodic_idx}.(branch).y(end), Manif.points{periodic_idx}.(branch).z(end)]; % last point
        P1 = [mappoints.x(idx-1), mappoints.y(idx-1), mappoints.z(idx-1)];
        P2 = [mappoints.x(idx),   mappoints.y(idx),   mappoints.z(idx)];
        first_return_dist = distance_to_line(P0, P1, P2);

        % fprintf('\n\n distance to line %e \n\n',first_return_dist);
        % Manif.grow_info.runinf.first_return_dist = first_return_dist;
% 
%         % check distance
%         if first_return_dist > abs(Manif.grow_info.init_step)
%             fprintf('Warning! the distance between last point of initial segment and its first return is greater than the initial step %e It is %e\n\n',abs(Manif.grow_info.init_step), first_return_dist)
%             prompt = "\n... Press something to continue\n\n";
%             x = input(prompt);
%         end
% 
        % make the new segment start at the end of the first intial segment
        mappoints.x = [Manif.points{periodic_idx}.(branch).x(end)  mappoints.x(idx:end)];
        mappoints.y = [Manif.points{periodic_idx}.(branch).y(end)  mappoints.y(idx:end)];
        mappoints.z = [Manif.points{periodic_idx}.(branch).z(end)  mappoints.z(idx:end)];
        idx_preimages(1:idx-1)=[];

        % %/////////////////
        % set(groot, 'CurrentFigure', fig_handles(idx_seg));
        % plot3(mappoints.x,mappoints.y,mappoints.z,'k.-');
        % %/////////////////
    end
% 
% 
%     if iter > period*mapiter  %if there is already a fdom on that branch
% 
%          %check the angles
%          % coordinates of mapped points
%          p0 = [Manif.points{idx_seg}.(branch).x(end-1), Manif.points{idx_seg}.(branch).y(end-1), Manif.points{idx_seg}.(branch).z(end-1)];
%          p1 = [mappoints.x(1), mappoints.y(1), mappoints.z(1)]; % the point we are actually looking at
%          p2 = [mappoints.x(2), mappoints.y(2), mappoints.z(2)];
%          delta0 = norm(p1-p0); % before
%          delta1 = norm(p1-p2); % after
%          alpha  = angles(p0,p1,p2); % angle btw points
% 
%         if (alpha*delta0 > manif.grow_info.deltalphamax || alpha*delta1 > manif.grow_info.deltalphamax )
%             fprintf('\n Warning! Delta*Alpha %e or %e between fund domains is larger than DeltaAlphaMax %e.\n Alpha=%e, Delta0=%e, Delta1=%e.\n The point is at x=%f y=%f z=%f\n\n', alpha*delta0, alpha*delta1, manif.grow_info.deltalphamax, alpha, delta0,delta1, mappoints.x(1), mappoints.y(1), mappoints.z(1)) 
%             prompt = "\n... Press something to continue\n\n";
%             x = input(prompt);
%         end
%     end
% 
%---%----------- Final section: save info

    %/////////////////
    % set(groot, 'CurrentFigure', fig_handles(periodic_idx));
    % plot3(mappoints.x,mappoints.y,mappoints.z,'.-','MarkerSize',15, 'LineWidth',1.5);
    %/////////////////

    mappoints.arc = arclength(mappoints);
    %chop the manifold up to the desired arclength (post mesh refinement)
    if mappoints.arc(end) > needed_arc
         idx_arc = find(mappoints.arc > needed_arc, 1, 'first'); %find first point that is more than the needed arc

         % !!! before used spline, now just chop at the point after!!!!
         % problems with arc duplicated up to numerical precision!!
         % chop mappoints up to the needed arclength
         mappoints.x   = mappoints.x(1:idx_arc);
         mappoints.y   = mappoints.y(1:idx_arc);
         mappoints.z   = mappoints.z(1:idx_arc);
         mappoints.arc = mappoints.arc(1:idx_arc);
         idx_preimages = idx_preimages(1:idx_arc);
    end


    

    % Add new segment to the entire manifold
    N = numel(Manif.points{periodic_idx}.(branch).x);
    Manif.points{periodic_idx}.(branch).x   = [Manif.points{periodic_idx}.(branch).x mappoints.x(2:end)]; %the first point of mappoints is the same as the last of previous fdom
    Manif.points{periodic_idx}.(branch).y   = [Manif.points{periodic_idx}.(branch).y mappoints.y(2:end)];
    Manif.points{periodic_idx}.(branch).z   = [Manif.points{periodic_idx}.(branch).z mappoints.z(2:end)];
    Manif.points{periodic_idx}.(branch).arc = arclength(Manif.points{periodic_idx}.(branch));

    Manif.points{periodic_idx}.(branch).idx_fdom(fdom_level+1,:) = [N numel(Manif.points{periodic_idx}.(branch).x)];
    Manif.points{periodic_idx}.(branch).idx_preimages     = [Manif.points{periodic_idx}.(branch).idx_preimages idx_preimages];

% 
% 
%     %Check once again if we reached the desired arclength
%     if strcmp(manif.orientability,'orientation-preserving') 
% 


    fprintf(' Fdomain arclength %.2f \n', mappoints.arc(end));
    fprintf(' %s %s branch arclength %.2f \n', Manif.points{periodic_idx}.name, branch, Manif.points{periodic_idx}.(branch).arc(end));

    if strcmp(branch,branch1)
        total_arc_branch1(periodic_idx) = Manif.points{periodic_idx}.(branch).arc(end);
    else
        total_arc_branch2(periodic_idx) = Manif.points{periodic_idx}.(branch).arc(end);
    end

    fprintf(' Total arclength of manifold: %.2f \n\n',  sum(total_arc_branch1)+sum(total_arc_branch2));


        
        % T{1,periodic_idx} = round(total_arc_branch1(periodic_idx), 2, 'significant');
        % disp(T)
%         %Check once again if we reached the desired arclength
%         if total_arc_branch1(idx_seg) > opts.user_arclength 
%             stop_arc_branch1(idx_seg) = 1;
%         end
% 
%     elseif strcmp(manif.orientability,'orientation-reversing')
%         if mod(floor(iter/period), 2) == 0 % 0: branch1, 1:branch2
% 
%             total_arc_branch1(idx_seg) = Manif.points{idx_seg}.(branch).arc(end);
%             fprintf(' Total arclength of manifold of %s %s: %.2f \n\n', Manif.points{idx_seg}.name, branch, total_arc_branch1(idx_seg));
%             T{1,idx_seg} = round(total_arc_branch1(idx_seg), 2, 'significant');
%             disp(T)
%             %Check once again if we reached the desired arclength
%             if total_arc_branch1(idx_seg) > opts.user_arclength 
%                 stop_arc_branch1(idx_seg) = 1;
%             end
% 
%         else
% 
%             total_arc_branch2(idx_seg) = Manif.points{idx_seg}.(branch).arc(end);
%             fprintf(' Total arclength of manifold of %s %s: %.2f \n\n',Manif.points{idx_seg}.name, branch, total_arc_branch2(idx_seg));
%             T{2,idx_seg} = round(total_arc_branch2(idx_seg), 2, 'significant');
%             disp(T)
%             %Check once again if we reached the desired arclength
%             if total_arc_branch2(idx_seg) > opts.user_arclength 
%                 stop_arc_branch2(idx_seg) = 1;
%             end
% 
%         end
%     end
% 
% % %---%----------- END Final section: save info
end
% % 
% % 
% % 
% % %--------------- END adding points
% % % 
% % % Erase last fundamental domain if the computation stopped chopping the last part of the manifold
% % if strcmp(manif.orientability,'orientation-preserving') && sum(stop_arc) == period
% %     for k=1:period
% %         Manif.points{k}.(branch1).idx_fund_dom(end, :) = [];
% %     end
% % elseif strcmp(manif.orientability,'orientation-reversing') && sum(stop_arc_branch1) + sum(stop_arc_branch2) == 2*period
% %     for k=1:period
% %         Manif.points{k}.(branch1).idx_fund_dom(end, :) = [];
% %         Manif.points{k}.(branch2).idx_fund_dom(end, :) = [];
% %     end
% % end
% % 
% % 
Manif.grow_info.runinf.time=toc;

fprintf('\n ----------------------')
if iter == manif.grow_info.max_iter
    fprintf('\n Algorithm stopping condition reached:\n Maximum number of iterations of Fdomains = %i ', iter);
elseif stop_arc == 1
    fprintf('\n Algorithm stopping condition reached:\n Arclength of the manifold has reached the (approx.) target arclength = %g ', opts.user_arclength);
end

fprintf('\n Elapsed time is %.3f seconds\n',Manif.grow_info.runinf.time)
fprintf(' ----------------------\n\n')
%%%%%%%%%



% progress table
varNames = arrayfun(@(i) sprintf('p%d_%d', period, i), order_periodic, 'UniformOutput', false);
if strcmp(manif.orientability,'orientation-preserving')
    data = num2cell(zeros(1, period), 1);
    T = table(data{:},'VariableNames',varNames,'RowName',{sprintf('%s arc%f\n',branch1)});
    T{1,:} = round(total_arc_branch1, 3, 'significant');
    disp(T)
    fprintf('    <strong>Total arc</strong> %g \n',  sum(total_arc_branch1));
else
    data = num2cell(zeros(2, period), 1);
    T = table(data{:},'VariableNames',varNames,'RowName',{sprintf('%s arc%f\n',branch1),sprintf('%s arc%f\n',branch2)});
    T{1,:} = round(total_arc_branch1, 2, 'significant');
    T{2,:} = round(total_arc_branch2, 2, 'significant');
    disp(T)
    fprintf('    <strong>Total arc</strong> %g \n',  sum(total_arc_branch1)+sum(total_arc_branch2));
end



fprintf('\n ----------------------')
fprintf('\n   * %i points removed \n',Manif.grow_info.runinf.rem_deltamin+Manif.grow_info.runinf.rem_nan+Manif.grow_info.runinf.rem_inf) 
fprintf('   * %i points added from deltamax \n',Manif.grow_info.runinf.add_deltamax) 
fprintf('   * %i points added from alpha \n',Manif.grow_info.runinf.add_alphamax) 
fprintf('   * %i points added from delta*alpha \n',Manif.grow_info.runinf.add_deltalphamax) 

fprintf(' ----------------------\n\n')

% % IF CONVERGENCE
% if sum(stop_convergence_branch1) == period 
%     fprintf('\n The arclength of all the last fdomains of the %s branch are less than the initial step |delta| = %.2e ', branch1, abs(manif.grow_info.init_step));
%     fprintf('\n The %s branch of the manifold appears to have reached convergence to: \n',branch1);
%     for k=1:period
%         fprintf('\n %s: (%.3f, %.3f, %.3f) ', Manif.points{k}.name, Manif.points{k}.(branch1).x(end), Manif.points{k}.(branch1).y(end), Manif.points{k}.(branch1).z(end));
%     end
%     fprintf('\n\n');
% elseif sum(stop_convergence_branch2) == period
%     fprintf('\n The arclength of all the last fdomains of the %s branch are less than the initial step |delta| = %.2e ', branch2, abs(manif.grow_info.init_step));
%     fprintf('\n The %s branch of the manifold appears to have reached convergence to: \n',branch2);
%     for k=1:period
%         fprintf('\n %s: (%.3f, %.3f, %.3f) ', Manif.points{k}.name, Manif.points{k}.(branch2).x(end), Manif.points{k}.(branch2).y(end), Manif.points{k}.(branch2).z(end));
%     end
%     fprintf('\n\n');
% end

%% FUNCTIONS

function interp = makima3D(points,t,tt)
% get interpolation points
% t: parametrization of the points
% tt: parametrization of the interpolated points
interp=struct();

interp.x = interp1(t,points.x,tt,'makima','extrap');
interp.y = interp1(t,points.y,tt,'makima','extrap');
interp.z = interp1(t,points.z,tt,'makima','extrap');
end 

%----------------

function arclen = arclength(points)
%arclength between each point of a vector (px,py,pz)
arclen=((points.x(1:end-1)-points.x(2:end)).^2 + (points.y(1:end-1)-points.y(2:end)).^2 + (points.z(1:end-1)-points.z(2:end)).^2).^(1/2);
arclen=[0 cumsum(arclen)];
end % function arclength

%----------------

function alpha=angles(p0,p1,p2)
%angle between p0p1 and p1p2
n1 = (p1 - p2) / norm(p1 - p2);  % Normalized vectors
n2 = (p0 - p1) / norm(p0 - p1);
alpha = atan2(norm(cross(n1, n2)), dot(n1, n2)); %gives value from 0 to pi
end

%----------------

function Anew=insert(A,B,ind)
% Anew: new vector with the new values
% A: Old vector
% B: vector with new values
%ind: index to insert values after this row


    % Preallocate output
    Anew = zeros(1,numel(A)+numel(B));

    % Find indices for old data
    addRows = ismember(1:numel(A), ind);
    oldDataInd = (1:numel(A)) + cumsum([0, addRows(1:end-1)]);

    % Add in old data
    Anew(oldDataInd) = A;

    % Find indices for new data
    newDataInd = (1:length(ind)) + ind;

    % Add in new data
    Anew(newDataInd) = B;
end

function dist = distance_to_line(P0, P1, P2)
    %distance Calculates the shortest distance from a point to a line in 3D.
    %   P0: A 1x3 vector representing the point (x0, y0, z0).
    %   P1: A 1x3 vector representing the first point on the line (x1, y1, z1).
    %   P2: A 1x3 vector representing the second point on the line (x2, y2, z2).
    %   dist: The shortest distance from P0 to the line defined by P1 and P2.
    
    
    % Vector from P1 to P0
    vec_P1P0 = P0 - P1;
    
    % Vector representing the line direction from P1 to P2
    vec_P1P2 = P2 - P1;
    
    % Calculate the cross product of vec_P1P0 and vec_P1P2
    cross_product = cross(vec_P1P0, vec_P1P2);
    
    % Calculate the magnitude (norm) of the cross product
    magnitude_cross_product = norm(cross_product);
    
    % Calculate the magnitude (norm) of the line direction vector
    magnitude_P1P2 = norm(vec_P1P2);

    % Calculate the distance using the formula
    dist = magnitude_cross_product / magnitude_P1P2;


end

end