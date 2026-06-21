function sol = solve_periodic_orbit(seed,opts,printt)

if nargin < 3
    printt = true;
end
    
    % --- 2. Guess matrix ---
    X0 = seed;
    n = size(X0, 1); % period


    % Set tolerances for high precision
    fsolve_options = optimoptions('fsolve', ...
    'Display', 'iter', ...
    'OptimalityTolerance', 1e-14, ... % Accuracy of the function value (residual)
    'StepTolerance', 1e-12, ...        % Accuracy of the change in X between steps
    'MaxFunctionEvaluations', 1e5, ... % Increase limits for complex n-period orbits
    'MaxIterations', 4000);

    if ~printt
        fsolve_options = optimoptions('fsolve','Display','none');
    end

    % --- 3. Solve the System ---
    % Target: find X such that P_{i+1} - f(P_i) = 0
        [sol, fval, exitflag] = fsolve(@(X) orbit_residual(X, n, opts), X0, fsolve_options);
    
    % --- 4. Process and Display Results ---
        if exitflag > 0 && printt
            orbit_points = sol;
            fprintf('\nFound a Period-%d Orbit:\n', n);
            disp(table((1:n)', orbit_points(:,1), orbit_points(:,2), orbit_points(:,3), ...
                'VariableNames', {'Orbit', 'x', 'y', 'z'}));
            fprintf('\nDistance to the original seed (Residuals):\n');
            diff = abs(orbit_points - X0);
    
            disp(table((1:n)', diff(:,1), diff(:,2), diff(:,3), ...
                'VariableNames', {'Orbit', 'dx', 'dy', 'dz'}));
        elseif exitflag <= 0 && printt
            fprintf('\nSolver failed to converge. Try a different seed.\n');
        end
    
    %%
    % --- The Residual Function (The 3n system) ---
    function F = orbit_residual(X, n, opts)
        thesystem=opts.thesystem;
        % Reshape vector back to a 3-by-n matrix for easy indexing
        F = zeros(n, 3);
        
        for i = 1:n
            % Define current point P_i
            current.x = X(i, 1);
            current.y = X(i, 2);
            current.z = X(i, 3);
            
            % Define the "target" point P_{i+1} 
            % (If i=n, the next point is 1, closing the loop)
            next_idx = mod(i, n) + 1;
            target = X(next_idx,:);
            
            % Map the current point forward once
            mapped = thesystem.ff(current, opts);
            
            % The equation: f(P_i) - P_{i+1} = 0
            F(i, 1) = mapped.x - target(1);
            F(i, 2) = mapped.y - target(2);
            F(i, 3) = mapped.z - target(3);
        end
    
    end
end