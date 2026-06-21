function orbit=eps_pseudo_orbit_periodic(manif, idxPO, idxpoint, branch)
% idxPO: which periodic point is the branch associated with (1 or 2 or ... period)
% idxpoint: the index to compute the pseudo orbit with
% branch: it can be branch1 or branch2
    
    i=0;

    while idxpoint~=0
        i=i+1;

        %where the point comes from
        orbit.name{i} = [manif.points{idxPO}.name branch];
        %coordinates of the point
        orbit.x(i) = manif.points{idxPO}.(branch).x(idxpoint);
        orbit.y(i) = manif.points{idxPO}.(branch).y(idxpoint);
        orbit.z(i) = manif.points{idxPO}.(branch).z(idxpoint);
        %index of the point in the manifold
        orbit.idxpoint(i) = idxpoint;

        % information for the next preimage
        idxpoint = manif.points{idxPO}.(branch).idx_preimages(idxpoint); %index of the manifold preimage
        branch   = manif.points{idxPO}.(branch).branch_preimage;         %branch of the preimage
        
        name     = manif.points{idxPO}.(branch).name_preimage;           %name of the preimage 
        idxPO    = str2double(name(strfind(name,'_')+1:end));            %index of periodic orbit of the preimage
        
    end


end