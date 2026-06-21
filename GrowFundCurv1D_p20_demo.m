%% Things to add:
% give an approximate and find with newton the true %% easiest
% ask for all periodic points of certain period %% a tiny bit harder

%% Adding the path
    clear all
    close all
    addpath('./GrowFundCurv1D_periodic_functions');

%% Options

    % %--- Define coordinates of the periodic orbit in original coordinates

       %orientation-preserving
       seed = ...
   [1.391783345971698   1.700271790018351   2.966957354816782;
   1.700271790018351   1.726610843859301   4.073837673871777;
   1.726610843859301   1.728896530872979   4.985680982956723;
   1.728896530872979   1.728900038693169   5.717441317238358;
   1.728900038693169   1.729573615468654   6.302853092483855;
   1.729573615468654   1.727245120282639   6.771856089455738;
   1.727245120282639   1.735496379100407   7.144729991847230;
   1.735496379100407   1.706225854214170   7.451280372578191;
   1.706225854214170   1.809442248141248   7.667250152276724;
   1.809442248141248   1.437786506905796   7.943242369962628;
   1.437786506905796   2.675602635002006   7.792380402875899;
   2.675602635002006  -2.527513508357938   8.909506957302725;
  -2.527513508357938  -1.385643744431251   4.600092057484241;
  -1.385643744431251   1.521737361011161   2.294429901556143;
   1.521737361011161   1.468622280773413   3.357281282256075;
   1.468622280773413   2.499669804719247   4.154447306578272;
   2.499669804719247  -1.607762448393132   5.823227649981866;
  -1.607762448393132   2.365000850952695   3.050819671592361;
   2.365000850952695  -1.875557759524913   4.805656588226583;
  -1.875557759524913   1.391783345971698   1.968967511056354];

  seed = ... %from sol 15000 
   [1.400640595070074   1.676256719002582   0.217685538728332;
   1.676256719002582   1.810355590519720   1.850405149985248;
   1.810355590519720   1.425489651574773   3.290679710507919;
   1.425489651574773   2.711085930409149   4.058033419981108;
   2.711085930409149  -2.722340026590012   5.957512666394035;
  -2.722340026590012  -2.397809441251363   2.043670106525217;
  -2.397809441251363  -2.366192124531177  -0.762873356031190;
  -2.366192124531177  -2.118208002568773  -2.976490809356129;
  -2.118208002568773  -0.996662779505743  -4.499400650053675;
  -0.996662779505743   2.571200903177254  -4.596183299548684;
   2.571200903177254  -2.710072918351250  -1.105745736461694;
  -2.710072918351250  -2.373134951827688  -3.594669507520605;
  -2.373134951827688  -2.244791375091578  -5.248870557844172;
  -2.244791375091578  -1.551028803233845  -6.443887821366917;
  -1.551028803233845   1.120872239011513  -6.706139060327379;
   1.120872239011513   2.478336782843165  -4.244039009250391;
   2.478336782843165  -1.605891537489955  -0.916894424557148;
  -1.605891537489955   2.364613404671097  -2.339407077135673;
   2.364613404671097  -1.873164014797224   0.493087742962558;
  -1.873164014797224   1.400640595070074  -1.478693820427177];
  % 

    %--- Information of the system
    opts.thesystem=StdHenon3D_periodic; % What is the name of the system file
    opts.par=struct('a',4.2,'b', -0.3, 'xi', 0.8); % The parameter values and names (has to match with the names defined in StdHenon3D)
    PO = solve_periodic_orbit(seed, opts);  %find periodic orbit using the seed and parallel shooting

    opts.user_arclength = 15000000; % What is the approximate arclength of the entire manifold
    opts.per_orbit.name ='p20';
    opts.per_orbit.coord = struct('x',PO(:,1)','y',PO(:,2)','z',PO(:,3)');
    opts.stability='Umanifold';

    %--- Number of iterations used to compute the manifold
    opts.max_iter =40*20;% 1000; % how many times (max) the algorithm iterates the fundamental domain
    opts.max_refines = 30;

    %--- Accuracy parameters (default)
    %opts.accpar.alphamax=0.3;
    %opts.accpar.deltalphamax=0.001; 
    %opts.accpar.deltamin=0.000001;
    %opts.accpar.deltamax=0.01;  

    %--- Initial step (default)
    opts.accpar.init_step=5*10^-3; % pos value is positive branch and viceversa

%% Computing the manifold: Wu(p7) orientation-preserving
     opts.branch = 'pos'; %which branch: 'pos', 'neg' or '' to consider sign of initial step.
     
     manif2 = GrowFundCurv1D_periodic(opts,2);

    % % adding the other branch of the manifold
     % manif = add_branch_periodic(manif, opts, 'neg');

%% Computing intersection points
    angle=pi/2; %the angle of the plane from [-pi, pi]. (angle=pi/2: x==0 (y>0), angle=0: y==0 (x>0))
    manif2=inter_plane_periodic(manif2,angle);


%% Plot
    % manifplot_periodic(manif2);

%% Epsilon pseudo orbit (orientation-preserving)
% 
% branch    = 'pos';
% idxPO     = 1;
% 
% if numel(manif.inter.points{idxPO}.(branch).idx)>0
% 
%     idxpoint  = manif.inter.points{idxPO}.(branch).idx(end);
%     orbit     = eps_pseudo_orbit_periodic(manif, idxPO, idxpoint, branch);
% 
%     % plot the epsilon orbit. Starting and end point are colored in solid red.
%     hold on
%     plot3(orbit.x([1,end]),orbit.y([1,end]),orbit.z([1,end]),'r.','MarkerSize',27) %epsilon orbit
%     plot3(orbit.x,orbit.y,orbit.z,'ko--','LineWidth',1.2) %epsilon orbit
% end
%%
% period = numel(PO(:,1));
% Ninter_pos=zeros(1,period);
% Ninter_neg=zeros(1,period);
% 
%     if isfield(manif2.points{1},'pos')
%         for k=1:period
%             Ninter_pos(k)=numel(manif2.inter.points{k}.pos.idx);
%         end
%     end
%     if isfield(manif2.points{1},'neg')
%         for k=1:period
%             Ninter_neg(k)=numel(manif2.inter.points{k}.neg.idx);
%         end
%     end
% 
% 
% Ninter_pos;
% Ninter_neg;
% Ntotal = sum(Ninter_pos) +  sum(Ninter_neg)
% log2(Ntotal)
% %%
% close all
% clear all
% 
% load("p20_15fdoms.mat")
% 
% period=numel(manif2.per_orbit.coord_original.x);
% 
% N_fdoms=[];
% if isfield(manif2.points{1},'pos')
%     for i=1:period
%         N_fdoms(end+1)= numel(manif2.points{i}.pos.idx_fdom(:,1));
%     end
% end
% 
% if isfield(manif2.points{1},'neg')
%     for i=1:period
%         N_fdoms(end+1)= numel(manif2.points{i}.neg.idx_fdom(:,1));
%     end
% end
% arc = zeros(1,period);
% N = min(N_fdoms);
% 
% 
% N=15;
% c = lines(period); %color
% 
% figure('Name', 'p20 x-coord');
% h1 = zeros(1, period);
% hold on
% 
% for per=1:period
%     h1(per) = subplot(period, 1, per);
%     hold(h1(per), 'on');
% end
% 
% if isfield(manif2.points{1},'pos')
%     branch = 'pos';
% 
%     for per = 1:period
%     c_idx = -per;
% 
%         for k=1:N%N_fdom
%             idx1 = manif2.points{per}.(branch).idx_fdom(k,1);
%             idx2 = manif2.points{per}.(branch).idx_fdom(k,2);
%             if manif2.points{per}.(branch).arc(idx2)-manif2.points{per}.(branch).arc(idx1)<1e-3
%                 plot(h1(per),manif2.points{per}.(branch).arc(idx1:idx2),manif2.points{per}.(branch).x(idx1:idx2),'*','color',c(mod(c_idx-1,period)+1, :));
%             else
%                 plot(h1(per),manif2.points{per}.(branch).arc(idx1:idx2),manif2.points{per}.(branch).x(idx1:idx2),'LineWidth',1.5,'color',c(mod(c_idx-1,period)+1, :));
%             end
%             c_idx = c_idx + 1;
%         end
%         arc(per) = max([arc(per) manif2.points{per}.(branch).arc(idx2)]);
%     end
% 
% end
% 
% if isfield(manif2.points{1},'neg')
%     branch = 'neg';
% 
%     for per = 1:period
%     c_idx = -per;
% 
%         for k=1:N%N_fdom
%             idx1 = manif2.points{per}.(branch).idx_fdom(k,1);
%             idx2 = manif2.points{per}.(branch).idx_fdom(k,2);
%             if manif2.points{per}.(branch).arc(idx2)-manif2.points{per}.(branch).arc(idx1)<1e-3
%                 plot(h1(per),-manif2.points{per}.(branch).arc(idx1:idx2),manif2.points{per}.(branch).x(idx1:idx2),'*','color',c(mod(c_idx-1,period)+1, :));
%             else
%                 plot(h1(per),-manif2.points{per}.(branch).arc(idx1:idx2),manif2.points{per}.(branch).x(idx1:idx2),'LineWidth',1.5,'color',c(mod(c_idx-1,period)+1, :));
%             end
%             c_idx = c_idx + 1;
%         end
%         arc(per) = max([arc(per) manif2.points{per}.(branch).arc(idx2)]);
%     end
% 
% end
% 
% 
% for per = 1:period
% 
%     h1(per) = subplot(period, 1, per);
%     plot(h1(per), 0, manif2.per_orbit.coord_compactified.x(per),'*','color',c(per, :))
% 
%     xline(0, 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1.5);
%     text(0, manif2.per_orbit.coord_compactified.x(per), manif2.per_orbit.name+"\_"+per);
%     ylabel(h1(per), sprintf('x'));
% 
%     ax = h1(per);
% 
%     if isfield(manif2.points{1},'pos') && isfield(manif2.points{1},'neg')
%         xlim(h1(per), [-arc(per) arc(per)]);
%         currentTicks = get(ax, 'XTick');
%         lastTick = max([-currentTicks(1) currentTicks(end)]);
%         deltaTick=currentTicks(2)-currentTicks(1);
%         % xlim(h1(per), [-lastTick-deltaTick lastTick+deltaTick]);
%         % currentTicks = get(ax, 'XTick');
%         % set(ax, 'XTick', sort(unique([currentTicks, -arc(per), 0, arc(per)])));
%     else
%         if isfield(manif2.points{1},'pos') 
%             xlim(h1(per), [0 arc(per)]);
%             currentTicks = get(ax, 'XTick');
%             lastTick = currentTicks(end);
%             deltaTick=currentTicks(2)-currentTicks(1);
%             xlim(h1(per), [0 lastTick+deltaTick]);
%         elseif isfield(manif2.points{1},'neg')   
%             xlim(h1(per), [-arc(per) 0]);
%             currentTicks = get(ax, 'XTick');
%             lastTick = currentTicks(1);
%             deltaTick=currentTicks(2)-currentTicks(1);
%             xlim(h1(per), [lastTick-deltaTick 0]);
%         end
%     end
% end
% 
% xlabel(h1(period), sprintf('arc'));
% 
% %%
% N=14;
% arcpos = zeros(1,period);
% 
% figure('Name', 'p20 x-coord pos');
% h2pos = zeros(1, period);
% hold on
% 
% for per=1:period
%     h2pos(per) = subplot(period, 1, per);
%     hold(h2pos(per), 'on');
% end
% 
% 
% branch = 'pos';
% 
% for per = 1:period
% c_idx = -per;
% 
%     for k=1:N%N_fdom
%         idx1 = manif2.points{per}.(branch).idx_fdom(k,1);
%         idx2 = manif2.points{per}.(branch).idx_fdom(k,2);
%         if manif2.points{per}.(branch).arc(idx2)-manif2.points{per}.(branch).arc(idx1)<1e-3
%             plot(h2pos(per),manif2.points{per}.(branch).arc(idx1:idx2),manif2.points{per}.(branch).x(idx1:idx2),'*','color',c(mod(c_idx-1,period)+1, :));
%         else
%             plot(h2pos(per),manif2.points{per}.(branch).arc(idx1:idx2),manif2.points{per}.(branch).x(idx1:idx2),'LineWidth',1.5,'color',c(mod(c_idx-1,period)+1, :));
%         end
%         c_idx = c_idx + 1;
%     end
%     arcpos(per) = max([arcpos(per) manif2.points{per}.(branch).arc(idx2)]);
% end
% 
% 
% 
% 
% for per = 1:period
% 
%     h2pos(per) = subplot(period, 1, per);
%     plot(h2pos(per), 0, manif2.per_orbit.coord_compactified.x(per),'*','color',c(per, :))
% 
%     xline(0, 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1.5);
%     text(0, manif2.per_orbit.coord_compactified.x(per), manif2.per_orbit.name+"\_"+per);
%     ylabel(h2pos(per), sprintf('x'));
% 
%     ax = h2pos(per);
% 
% 
%     xlim(h2pos(per), [0 arcpos(per)]);
%     currentTicks = get(ax, 'XTick');
%     lastTick = currentTicks(end);
%     deltaTick=currentTicks(2)-currentTicks(1);
%     % xlim(h2pos(per), [0 lastTick+deltaTick]);
% end
% xlabel(h2pos(period), sprintf('arc'));
% 
% %%
% arcneg = zeros(1,period);
% 
% figure('Name', 'p20 x-coord neg');
% h2neg = zeros(1, period);
% hold on
% 
% for per=1:period
%     h2neg(per) = subplot(period, 1, per);
%     hold(h2neg(per), 'on');
% end
% 
% branch = 'neg';
% 
% for per = 1:period
% c_idx = -per;
% 
%     for k=1:N%N_fdom
%         idx1 = manif2.points{per}.(branch).idx_fdom(k,1);
%         idx2 = manif2.points{per}.(branch).idx_fdom(k,2);
%         if manif2.points{per}.(branch).arc(idx2)-manif2.points{per}.(branch).arc(idx1)<1e-3
%             plot(h2neg(per),-manif2.points{per}.(branch).arc(idx1:idx2),manif2.points{per}.(branch).x(idx1:idx2),'*','color',c(mod(c_idx-1,period)+1, :));
%         else
%             plot(h2neg(per),-manif2.points{per}.(branch).arc(idx1:idx2),manif2.points{per}.(branch).x(idx1:idx2),'LineWidth',1.5,'color',c(mod(c_idx-1,period)+1, :));
%         end
%         c_idx = c_idx + 1;
%     end
%     arcneg(per) = max([arcneg(per) manif2.points{per}.(branch).arc(idx2)]);
% end
% 
% 
% 
% 
% for per = 1:period
% 
%     h2neg(per) = subplot(period, 1, per);
%     plot(h2neg(per), 0, manif2.per_orbit.coord_compactified.x(per),'*','color',c(per, :))
% 
%     xline(0, 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1.5);
%     text(0, manif2.per_orbit.coord_compactified.x(per), manif2.per_orbit.name+"\_"+per);
%     ylabel(h2neg(per), sprintf('x'));
% 
%     ax = h2neg(per);
% 
% 
%     xlim(h2neg(per), [-arcneg(per) 0]);
%     currentTicks = get(ax, 'XTick');
%     lastTick = currentTicks(end);
%     deltaTick=currentTicks(2)-currentTicks(1);
%     % xlim(h2pos(per), [0 lastTick+deltaTick]);
% end
% xlabel(h2neg(period), sprintf('arc'));
% 
% 
% %% plot 3D
% per_orbit = manif2.per_orbit.coord_compactified;
% green=[0 204 26]/255;
% mks=11;
% lnw=1.7;
% 
% h = figure;
% hold on
% box on
% grid on
% 
% for per=1:period
%     c_idx = -per;
%     for k=1:N%N_fdom
%         idx1 = manif2.points{per}.(branch).idx_fdom(k,1);
%         idx2 = manif2.points{per}.(branch).idx_fdom(k,2);
%         plot3(manif2.points{per}.(branch).x(idx1:idx2),manif2.points{per}.(branch).y(idx1:idx2),manif2.points{per}.(branch).z(idx1:idx2),'LineWidth',lnw,'color',c(mod(c_idx-1,period)+1, :));
%         c_idx = c_idx + 1;
%     end
% end
% view([290 30])
% 
% plot3(per_orbit.x,per_orbit.y,per_orbit.z,'.','marker','o','MarkerFaceColor',green,'MarkerEdgeColor',green,'LineWidth',1,'MarkerSize',5)
% 
% %-- Unit circle
% [xunit,yunit] = circle(0,0,1,1000);
% plot3(xunit,yunit,ones(size(xunit)),'k','LineWidth',1.5)
% plot3(xunit,yunit,-ones(size(xunit)),'k','LineWidth',1.5)
% % % 
% % view([290 30])
% % % daspect([1 1 1])
% % allx=[manif2.points{1}.(branch).x manif2.points{2}.(branch).x manif2.points{3}.(branch).x];
% % ally=[manif2.points{1}.(branch).y manif2.points{2}.(branch).y manif2.points{3}.(branch).y];
% % allz=[manif2.points{1}.(branch).z manif2.points{2}.(branch).z manif2.points{3}.(branch).z];
% % 
% % ax = gca;
% % ax.XTick = [-0.5 0 0.5];
% % ax.YTick = [-0.5 0 0.5];
% % ax.ZTick = [0.6 0.7 0.8];
% % 
% % xlim([min(allx)-.1 max(allx)])
% % ylim([min(ally) max(ally)])
% % zlim([min(allz) max(allz)])
% 

%%



% > --------  eigenvalue and eigenvector in compactified coordinates
function [eigvecs_comp, eigval] = eigensystem(per_orbit, opts, inv_flag)
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
    syms x y z
    vars = [x, y, z];
    points = struct('x',x, 'y',y, 'z',z);
    system = opts.thesystem;
    period = numel(per_orbit.x);
    n = 3; % Dim   

    %Jacobian J_f of the original system (uncompactified)
    if inv_flag == 1
        F = system.ff_inv(points,opts);
    else
        F = system.ff(points,opts);
    end
    %Jacobian J_t of the compactification
    T = system.compactify(points);


    % Convert Symbolic to Fast Numeric Functions
    JF_func = matlabFunction(jacobian([F.x, F.y, F.z], vars), 'Vars', {vars});
    JT_func = matlabFunction(jacobian([T.x, T.y, T.z], vars), 'Vars', {vars});


    %Pre-compute all local Jacobians and the Jacobian of the k-th iterate
    %Jacobian of the k-th iterate: J_f^k(xi)= J_f(xi-1) * J_f(xi-2))* ... * J_f(x0) * J_f(xk) * ... * J_f(xi+1) * J_f(xi) 
    % Storing them in a 3D array is memory-efficient for small n
    JF_all = zeros(n, n, period); % initialize. all_JF(:,:i) is the local Jacobian of xi
    JFk = eye(n); %identity
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
    eigvecs_comp  = zeros(n, n, period);

    %%%%
    for i = 1:period
            pt = [per_orbit.x(i), per_orbit.y(i), per_orbit.z(i)];
            
            % Evaluate functions
            JTi = JT_func(pt);

            % Transform current vectors to compact space and normalize
            V_comp = JTi * VF_current; %compactify last computed eigenvector
            eigvecs_comp(:,:,i) = V_comp ./ vecnorm(V_comp, 2, 1); %normalize
    
            % Propagate original eigenvectors (original coord) to the NEXT point in the orbit
            if i < period
                VF_next = JF_all(:,:,i) * VF_current; %v_{i+1} = Jf_i * v_i
                % update:
                VF_current = VF_next ./ vecnorm(VF_next, 2, 1); %normalize
            end
    end

end

%----------------

function arclen = arclength(points)
%arclength between each point of a vector (px,py,pz)
arclen=((points.x(1:end-1)-points.x(2:end)).^2 + (points.y(1:end-1)-points.y(2:end)).^2 + (points.z(1:end-1)-points.z(2:end)).^2).^(1/2);
arclen=[0 cumsum(arclen)];
end % function arclength

%----------------
