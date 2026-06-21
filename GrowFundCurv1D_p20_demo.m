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

  %orientation-reversing
  seed = ...
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

    opts.user_arclength = 100000; % What is the approximate arclength of the entire manifold
    opts.per_orbit.name ='p20';
    opts.per_orbit.coord = struct('x',PO(:,1)','y',PO(:,2)','z',PO(:,3)');
    opts.stability='Umanifold';

    %--- Number of iterations used to compute the manifold
    opts.max_iter =40*10; % how many times (max) the algorithm iterates the fundamental subdomains
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
     
     manif = GrowFundCurv1D_periodic(opts,2);

    % % adding the other branch of the manifold
     % manif = add_branch_periodic(manif, opts, 'neg');

%% Computing intersection points
    angle=pi/2; %the angle of the plane from [-pi, pi]. (angle=pi/2: x==0 (y>0), angle=0: y==0 (x>0))
    manif=inter_plane_periodic(manif,angle);


%% Plot
    % manifplot_periodic(manif);

%----------------
