function plotState(surfs)
        nSurfs = length(surfs);
        %ColorSet = varycolor(Nb*(N_span+1));
        figure(101); clf; hold all;
        %set(gca, 'ColorOrder', ColorSet);
%         for n = 1:Nb
%             quiver3(x_cp(n,:),y_cp(n,:),z_cp(n,:),...
%                     control_normals(n,:,1),...
%                     control_normals(n,:,2),...
%                     control_normals(n,:,3));
%         end

        for n = 1:nSurfs
           sSurf = size(surfs(n).xSurface);
           sWake = size(surfs(n).xWake);
           N_span = sSurf(1);
           N_chord = sSurf(2);
           N_wake = sWake(2);
           for i = 1:N_span
              plot3(surfs(n).xSurface(i,:),...
                    surfs(n).ySurface(i,:),...
                    -surfs(n).zSurface(i,:),'-k')
           end
           for j = 1:N_chord
              plot3(surfs(n).xSurface(:,j),...
                    surfs(n).ySurface(:,j),...
                    -surfs(n).zSurface(:,j),'-k')
           end
           for i = 1:N_span
              plot3(surfs(n).xWake(i,:),...
                    surfs(n).yWake(i,:),...
                    -surfs(n).zWake(i,:), '-b')
           end
           for j = 1:N_wake
              plot3(surfs(n).xWake(:,j),...
                    surfs(n).yWake(:,j),...
                    -surfs(n).zWake(:,j), '-b')
           end
        end
        axis equal;
        view([0,1,0])
    end