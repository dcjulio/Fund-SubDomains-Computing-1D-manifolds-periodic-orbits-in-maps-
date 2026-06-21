function h=manifplot_periodic(manif)


    per_orbit = manif.per_orbit.coord_compactified;
    period    = numel(manif.points);

    if strcmp(manif.stability,'Smanifold')
        RGB1=[131, 195, 251]/255; % light blue
        RGB2=[5, 52, 122]/255; % dark blue
    elseif strcmp(manif.stability,'Umanifold')
        RGB1=[240, 120, 98]/255; % light red
        RGB2=[168, 25, 17]/255;
    end
    R= linspace(RGB1(1),RGB2(1),period);  %// Red from 212/255 to 0
    G = linspace(RGB1(2),RGB2(2),period);   %// Green from 212/255 to 0
    B = linspace(RGB1(3),RGB2(3),period);  %// Blue from 1 to 170/255
    c = [R(:), G(:), B(:)];



    h=figure;
    plot3(per_orbit.x,per_orbit.y,per_orbit.z,'*')
    hold on
    
%%
   %plot each branch
    if isfield(manif.points{1},'pos')
        for i=1:period
            h = plot3(manif.points{i}.pos.x,manif.points{i}.pos.y,manif.points{i}.pos.z,'LineWidth',2,'color',c(i,:));
            h.UserData.branch_name = strrep([manif.points{i}.name 'pos'], '_', '\_');
        end

        % Enable data cursor mode
        dcm = datacursormode(gcf);
        dcm.Enable = 'on';
        % Custom update function for the data tip
        dcm.UpdateFcn = @(obj, event) myCustomDataTip(obj, event);
    end


    if isfield(manif.points{1},'neg')
        for i=1:period
            h = plot3(manif.points{i}.neg.x,manif.points{i}.neg.y,manif.points{i}.neg.z,'LineWidth',2,'color',c(i,:));
            h.UserData.branch_name = strrep([manif.points{i}.name 'neg'], '_', '\_');
        end
        % Enable data cursor mode
        dcm = datacursormode(gcf);
        dcm.Enable = 'on';
        % Custom update function for the data tip
        dcm.UpdateFcn = @(obj, event) myCustomDataTip(obj, event);
    end

%%
    % % To plot each fundamental domain in different colour
    % if isfield(manif.points{1},'pos')
    %     for i=1:period
    %         nRows = size(manif.points{i}.pos.idx_fund_dom, 1);
    %         idxs=[0 0];
    %         for r=1:nRows
    %             idxs=manif.points{i}.pos.idx_fund_dom(r,:);
    %             h=plot3(manif.points{i}.pos.x(idxs(1):idxs(2)),manif.points{i}.pos.y(idxs(1):idxs(2)),manif.points{i}.pos.z(idxs(1):idxs(2)),'LineWidth',2);
    %             branch_name=sprintf('%s pos fund_dom: %i ',manif.points{i}.name,r);
    %             h.UserData.branch_name = strrep(branch_name, '_', '\_');
    %         end
    %             h=plot3(manif.points{i}.pos.x(idxs(2)+1:end),manif.points{i}.pos.y(idxs(2)+1:end),manif.points{i}.pos.z(idxs(2)+1:end),'LineWidth',2);
    %             branch_name=sprintf('%s pos fund_dom: %i ',manif.points{i}.name,r+1);
    %             h.UserData.branch_name = strrep(branch_name, '_', '\_');
    %     end
    %     % Enable data cursor mode
    %     dcm = datacursormode(gcf);
    %     dcm.Enable = 'on';
    %     % Custom update function for the data tip
    %     dcm.UpdateFcn = @(obj, event) myCustomDataTip(obj, event);
    % 
    % end
    % if isfield(manif.points{1},'neg')
    %     for i=1:period
    %         nRows = size(manif.points{i}.neg.idx_fund_dom, 1);
    %         idxs=[0 0];
    %         for r=1:nRows
    %             idxs=manif.points{i}.neg.idx_fund_dom(r,:);
    %             h=plot3(manif.points{i}.neg.x(idxs(1):idxs(2)),manif.points{i}.neg.y(idxs(1):idxs(2)),manif.points{i}.neg.z(idxs(1):idxs(2)),'LineWidth',2);
    %             branch_name=sprintf('%s neg fund_dom: %i ',manif.points{i}.name,r);
    %             h.UserData.branch_name = strrep(branch_name, '_', '\_');
    %         end
    %             h=plot3(manif.points{i}.neg.x(idxs(2)+1:end),manif.points{i}.neg.y(idxs(2)+1:end),manif.points{i}.neg.z(idxs(2)+1:end),'LineWidth',2);
    %             branch_name=sprintf('%s neg fund_dom: %i ',manif.points{i}.name,r+1);
    %             h.UserData.branch_name = strrep(branch_name, '_', '\_');
    %     end
    %     % Enable data cursor mode
    %     dcm = datacursormode(gcf);
    %     dcm.Enable = 'on';
    %     % Custom update function for the data tip
    %     dcm.UpdateFcn = @(obj, event) myCustomDataTip(obj, event);
    % end

    %%
    for k=1:period
        text(per_orbit.x(k),per_orbit.y(k),per_orbit.z(k),sprintf('%s', strrep(manif.points{k}.name,'_','\_')))
    end

    
    if isfield(manif,'inter')
        angle=str2num(manif.inter.angle(1:end-2))*pi;
        % % %-- intersection points
        if isfield(manif.points{1},'pos')

            for k = 1:period
                x=manif.inter.points{k}.pos.x;
                y=manif.inter.points{k}.pos.y;
                z=manif.inter.points{k}.pos.z;

                if numel(x)>0
                    h = plot3(x,y,z,'.','color',RGB2*0.7,'MarkerSize',11);
                    % Attach custom data
                    h.UserData.branch_name = strrep([manif.points{k}.name 'pos'], '_', '\_');
                end

            end
            
            % Enable data cursor mode
            dcm = datacursormode(gcf);
            dcm.Enable = 'on';
            % Custom update function for the data tip
            dcm.UpdateFcn = @(obj, event) myCustomDataTip(obj, event);

        end

        if isfield(manif.points{1},'neg')
            for k = 1:period
                x=manif.inter.points{k}.neg.x;
                y=manif.inter.points{k}.neg.y;
                z=manif.inter.points{k}.neg.z;

                if numel(x)>0
                    h = plot3(x,y,z,'.','color',RGB2*0.7,'MarkerSize',11);
                    % Attach custom data
                    h.UserData.branch_name = strrep([manif.points{k}.name 'neg'], '_', '\_');
                end
            end
            
            % Enable data cursor mode
            dcm = datacursormode(gcf);
            dcm.Enable = 'on';
            % Custom update function for the data tip
            dcm.UpdateFcn = @(obj, event) myCustomDataTip(obj, event);

        end

        %-- Plane
        plane.x=[cos(angle),0];
        plane.y=[sin(angle),0];
        plane.z= [-1,1];
        plane.color=[230, 178, 17]/255;
        surf(repmat(plane.x,2,1), repmat(plane.y,2,1), repmat(plane.z,2,1)','FaceAlpha',0.4, 'EdgeColor',plane.color,'FaceColor',plane.color,'FaceLighting','gouraud','LineWidth',1.7)

    end


%-- Unit circle
[xunit,yunit] = circle(0,0,1,1000);
plot3(xunit,yunit,ones(size(xunit)),'k','LineWidth',1.5)
plot3(xunit,yunit,-ones(size(xunit)),'k','LineWidth',1.5)
% 
       
xlabel('x')
ylabel('y')
zlabel('z')
xlim([-1.01 1.01])
ylim([-1.01 1.01])
zlim([-1.01 1.01])


title(strrep(manif.name,'_','\_'))

daspect([1 1 1])
view([100,30])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xunit,yunit] = circle(x,y,r,n)
th = linspace(0,2*pi,n);
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
end

function txt = myCustomDataTip(~, event_obj)
    pos = event_obj.Position;
    idx = event_obj.DataIndex;
    target = event_obj.Target;
    
    % Get the k value from UserData
    if isfield(target.UserData, 'branch_name')
        kk = target.UserData.branch_name;
    else
        kk = NaN;
    end
    

    % Custom text
    txt = {
        ['X: ', num2str(pos(1))], ...
        ['Y: ', num2str(pos(2))], ...
        ['Z: ', num2str(pos(2))], ...
        ['idx: ', num2str(idx)], ...
        ['Branch name: ', kk]
    };
end

end