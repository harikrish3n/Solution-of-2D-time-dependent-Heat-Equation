% Class containing additional utilities such as displating multiple surface
% plots and generating the plots in the required format
classdef Utilities
    methods(Static)

        function surface_plot_expl(T,Nx,Ny,dt,time)
            hx = 1./(Nx+1); hy = 1./(Ny+1);
            X = cell(size(Nx,2)); Y = cell(size(Ny,2));
            for i=1:size(T,1)
                [X{i},Y{i}] = meshgrid(0:hx(i):1,0:hy(i):1);
                for j=1:size(T,2)
                    subplot(size(T,1),size(T,2),j + (i-1)*size(dt,2));
                    surf(X{i},Y{i},T{i,j},'EdgeColor','interp');
                    title(['Nx,Ny = ' num2str(Nx(i)) '; dt = 1/' num2str(1/dt(j))]);
                    xlabel("X");ylabel('Y'); zlabel('T');
                end
            end
            sgtitle(['Solutions using explicit Euler Method at time = ' num2str(time*8) '/8']);
        end

        function surface_plot_impl(T,Nx,Ny,dt,time)
            hx = 1./(Nx+1); hy = 1./(Ny+1);
            X = cell(size(Nx,2)); Y = cell(size(Ny,2));
            for i=1:size(T,1)
                [X{i},Y{i}] = meshgrid(0:hx(i):1,0:hy(i):1);
                subplot(1,4,i);
                surf(X{i},Y{i},T{i},'EdgeColor','interp');
                title(['Nx,Ny = ' num2str(Nx(i)) '; dt = 1/' num2str(1/dt)]);
                xlabel("X");ylabel('Y'); zlabel('T');
            end
            sgtitle(['Solutions using implicit Euler Method at time = ' num2str(time*8) '/8']);
        end

        function generate_plot(T,Nx,Ny,dt,time)
            hx = 1./(Nx+1); hy = 1./(Ny+1);
            X = cell(size(Nx,2)); Y = cell(size(Ny,2));
            for i=1:size(T,1)
                [X{i},Y{i}] = meshgrid(0:hx(i):1,0:hy(i):1);
                for j=1:size(T,2)
                    f = figure('Visible','off');
                    surf(X{i},Y{i},T{i,j},'EdgeColor','interp');
                    title(['Explicit Euler Method at time =' num2str(time*8) '/8'],['Nx,Ny = ' num2str(Nx(i)) '; dt = 1/' num2str(1/dt(j))]);
                    xlabel("X");ylabel('Y'); zlabel('T');
                    exportgraphics(f,"NXY"+Nx(i)+"_DT1_"+1/dt(j)+"_TIME"+time*8+"_8.png",'Resolution',300)
                end
            end
        end
  
    end
end