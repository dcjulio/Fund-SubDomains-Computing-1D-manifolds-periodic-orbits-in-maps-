function manif = init_manif_periodic(opts)

%---- manif.name: Name of the manifold---% 
% for example Ws_pmin_pos is the stable manifolds of the fixed point pmin, the branch to the positive values
%
%---- manif.orientability: Orientability of the manifold---% 
%
%---- manif.fixp: Information of the fixed points associated to the manifold---% 
%
%---- manif.stability: Stability of the manifold---% 
%
%---- manif.points: Coordinates of the manifold ---% 
%
%---- manif.system_info: contains general info about the map ---% 
% system_info.par: The parameter values
% system_info.fixp: the fixed points with their  eigensystem, orientability, stability, etc
%
%---- manif.grow_info: information for the algorithm---% 
%
%% to-do
% add syst_info with fixed points if defined in the class of functions. % eigensystem, stability and orientability of all the fixed points

%%
%the map function where the system is defined (example StdHenon3D)
thesystem=opts.thesystem;   
%% Initializate field names and the structure 'manif'

%names of the fields
names = {
    'name'
    'par' % parameters
    'type' %periodic orbit or fixed point
    'orientability'
    'stability'   % stability
    'dimension'
    'per_orbit'   % periodic orbit
    'grow_info'    % options for computation
    }; 

manif=struct();
n = numel(names);

for k = 1:n
    manif.(names{k})=[];
end

manif.type = sprintf('period-%i orbit', numel(opts.per_orbit.coord.x));
%% Default accuracty conditions

% default acc conditions
manif.grow_info.alphamax=0.3;
manif.grow_info.deltalphamax=10^(-3); 
manif.grow_info.deltamin=10^(-6);
manif.grow_info.deltamax=10^(-2);    
manif.grow_info.init_step=10^(-7);  

%% Update accuracy conditions if needed
%rewrite them if we other values are defined
manif.grow_info.thesystem = thesystem;

if isfield(opts,'accpar')
    if isfield(opts.accpar,'alphamax')
        manif.grow_info.alphamax=opts.accpar.alphamax;
    end
    if isfield(opts.accpar,'deltalphamax')
        manif.grow_info.deltalphamax=opts.accpar.deltalphamax;
    end
    if isfield(opts.accpar,'deltamin')
        manif.grow_info.deltamin=opts.accpar.deltamin;
    end
    if isfield(opts.accpar,'deltamax')
        manif.grow_info.deltamax=opts.accpar.deltamax;
    end
    if isfield(opts.accpar,'init_step')
        manif.grow_info.init_step=opts.accpar.init_step;
    end
end

%% Parameter values
manif.par=opts.par;  % parameters

%% General info of periodic orbit
manif.per_orbit.name=opts.per_orbit.name;
manif.per_orbit.coord_original=opts.per_orbit.coord;
manif.per_orbit.coord_compactified=thesystem.compactify(opts.per_orbit.coord); %periodic orbit in compactified coordinates

% eigensystem
[manif.per_orbit.eigvecs_comp, manif.per_orbit.eigval, manif.per_orbit.mu] = eigensystem(opts.per_orbit.coord,opts,0);
eigval=manif.per_orbit.eigval;
eigvecs=manif.per_orbit.eigvecs_comp;

% computes the dimension and orientation properties of the stable manifold
manif.per_orbit.Smanifold.dimension=sum(abs(eigval)<1);
manif.per_orbit.Smanifold.orientability=orientability(eigval,'Smanifold',opts);
manif.per_orbit.Smanifold.eigval=eigval(abs(eigval)<1);
manif.per_orbit.Smanifold.eigvecs=eigvecs(:,abs(eigval)<1,:);
if strcmp(opts.stability,'Smanifold')
    manif.dimension=manif.per_orbit.Smanifold.dimension;
    manif.orientability=manif.per_orbit.Smanifold.orientability;
end

% computes the dimension and orientation properties of the unstable manifold
manif.per_orbit.Umanifold.dimension=sum(abs(eigval)>1);
manif.per_orbit.Umanifold.orientability=orientability(eigval,'Umanifold',opts);
manif.per_orbit.Umanifold.eigval=eigval(abs(eigval)>1);
manif.per_orbit.Umanifold.eigvecs=eigvecs(:,abs(eigval)>1,:);
if strcmp(opts.stability,'Umanifold')
    manif.dimension=manif.per_orbit.Umanifold.dimension;
    manif.orientability=manif.per_orbit.Umanifold.orientability;
end


%% Stability and orientability and dimension of the manif to compute

manif.stability=opts.stability;

%if the field branch is defined, then follow that definition to know which branch to compute
if isfield(opts,'branch')
   if strcmp(opts.branch,'pos')
       manif.grow_info.init_step=abs(manif.grow_info.init_step);
   elseif strcmp(opts.branch,'neg')
       manif.grow_info.init_step=-abs(manif.grow_info.init_step);
   end
end

%% Name of the manifold. Example: Ws_pmin

% defining the name of the manifold
manif.name = sprintf('W%s_%s', lower(manif.stability(1)),opts.per_orbit.name);  

%% algorithm information
manif.grow_info.stability=manif.stability; % stability of the manifold
manif.grow_info.orientability=manif.orientability; % orientability of the manifold
manif.grow_info.dimension=manif.dimension; % dimension of the manifold

manif.grow_info.eigval=manif.per_orbit.(manif.stability).eigval; 
manif.grow_info.eigvecs=manif.per_orbit.(manif.stability).eigvecs; 

manif.grow_info.max_iter=opts.max_iter; % max number of iteration of fdoms
manif.grow_info.user_arclength=opts.user_arclength;

   
%----------------------------------------------
%-------------- FUNCTIONS ---------------------
%----------------------------------------------

% > --------  eigenvalue and eigenvector in compactified coordinates
function [eigvecs_comp, eigval, mu] = eigensystem(per_orbit, opts, inv_flag)
% --- Inputs ---
    % per_orbit: struct with fields .x, .y, .z (vectors of length period = k)
    % opts: system options
    % inv_flag: 1 for inverse, 0 for forward
% --- Outputs ---
    % eigval: the eigenvalue of the periodic orbit
    % eigvecs(:,k,period): each column k is an eigenvector associated to the kth eigenval
% ---
% ---
    % Creating Symbolic Variables and Jacobians
    vars_sym = sym('v', [1 3]); % Create a symbolic vector
    x_sym = vars_sym(1); y_sym = vars_sym(2); z_sym = vars_sym(3);

    points = struct('x', x_sym, 'y', y_sym, 'z', z_sym);
    system = opts.thesystem;
    period = numel(per_orbit.x);
    dim = 3; % Dim   

    %Jacobian J_f of the original system (uncompactified)
    if inv_flag == 1
        F = system.ff_inv(points,opts);
    else
        F = system.ff(points,opts);
    end
    %Jacobian J_t of the compactification
    T = system.compactify(points);


    % Convert Symbolic to Fast Numeric Functions
    JF_func = matlabFunction(jacobian([F.x, F.y, F.z], vars_sym), 'Vars', {vars_sym});
    JT_func = matlabFunction(jacobian([T.x, T.y, T.z], vars_sym), 'Vars', {vars_sym});


    %Pre-compute all local Jacobians and the Jacobian of the k-th iterate
    %Jacobian of the k-th iterate: J_f^k(xi)= J_f(xi-1) * J_f(xi-2))* ... * J_f(x0) * J_f(xk) * ... * J_f(xi+1) * J_f(xi) 
    % Storing them in a 3D array is memory-efficient for small n
    JF_all = zeros(dim, dim, period); % initialize. all_JF(:,:i) is the local Jacobian of xi
    JFk = eye(dim); %identity
    for i = 1:period
        pt = [per_orbit.x(i), per_orbit.y(i), per_orbit.z(i)];
        JFi = JF_func(pt);
        JF_all(:,:,i) = JFi;
        JFk = JFi * JFk; % Cumulative product (Jacobian of the k-th iterate at x0)
    end

    % Stability of the Orbit (at x0)
    [VF_current, D] = eig(JFk);
    eigval = diag(D).'; % eigenvalues of JF^k
    %eigvals = zeros(period, n); % To store local eigenvalues if desired

    % Propagate Eigenvectors and Transform to Compactified Coordinates
    % We store all eigenvectors for all points in a 3D array: [dim x num_vectors x period]
    eigvecs_comp  = zeros(dim, dim, period);


    %%%%
    for i = 1:period
            pt = [per_orbit.x(i), per_orbit.y(i), per_orbit.z(i)];
            
            % Evaluate functions
            JTi = JT_func(pt);

            % Transform current vectors to compact space and normalize
            V_comp = JTi * VF_current; %compactify last computed eigenvector

            eigvecs_comp(:,:,i) = V_comp ./ vecnorm(V_comp, 2, 1); %normalize
    
             VF_next = JF_all(:,:,i) * VF_current; %v_{i+1}*mui = Jf_i * v_i
             % update:
             mu(:,i)=vecnorm(VF_next, 2, 1);

            % Propagate original eigenvectors (original coord) to the NEXT point in the orbit
            if i < period
                VF_current = VF_next ./ mu(:,i)'; %normalize
            end

    end

end
%----------------------------------------------
%----------------------------------------------
%----------------------------------------------


% > -------- orientability
function [orientability]=orientability(eigval,Stab,opts)
    
    % If stable manifold, the eigenvalue is less than one and its more numerically unstable. 
    % Hence, we use the inverserse to get the orientability, less numerical innacuracy to get the sign
    if strcmp(Stab,'Smanifold') 
            inv=1; % use the inverse
            per_orbit = opts.per_orbit.coord;
            [~,eigval,~]=eigensystem(per_orbit,opts,inv);

            if prod(eigval(abs(eigval) > 1)) > 0
                orientability='orientation-preserving';
            else
                orientability='orientation-reversing';
            end
    end


    if strcmp(Stab,'Umanifold')
        if prod(eigval(abs(eigval) > 1)) > 0
            orientability='orientation-preserving';
        else
            orientability='orientation-reversing';
        end
    end
    
end  

% > -------- oarclength
function arclen = arclength(points)
%arclength between each point of a vector (px,py,pz)

arclen=((points.x(1:end-1)-points.x(2:end)).^2 + (points.y(1:end-1)-points.y(2:end)).^2 + (points.z(1:end-1)-points.z(2:end)).^2).^(1/2);
arclen=[0 cumsum(arclen)];
end % function arclength
    
end