% Class containing the implementation of the explicit and implicit Euler
% methods along with the stability check

classdef Numerical_Methods
    methods(Static)

        function T_new = expl_euler(Nx,Ny,dt,T)
            hx=1/(Nx+1);hy=1/(Ny+1);
            T_new = T; % Can also initialize with zeros since it is explicit
            for j=2:Ny+1
                for i=2:Nx+1
                    T_new(i,j) =T(i,j) + dt*(((T(i-1,j) - 2*T(i,j) + T(i+1,j))/hx^2) + ((T(i,j-1) - 2*T(i,j) + T(i,j+1))/hy^2));
                end
            end
        end

        function T_new = impl_euler(Nx,Ny,dt,T)
            hx=1/(Nx+1);hy=1/(Ny+1);
            T_new = T;
            norm=1;
            while(norm>1e-6)
                norm=0;
                for j=2:Ny+1
                    for i=2:Nx+1
                        T_new(i,j) = (T(i,j) + (dt/hx^2)*(T_new(i-1,j) + T_new(i+1,j)) + (dt/hy^2)*(T_new(i,j-1) + T_new(i,j+1)))/(1 + (2*dt/hx^2) + (2*dt/hy^2));
                    end
                end
                for j=2:Ny+1
                    for i=2:Nx+1
                        norm = norm + ((T_new(i,j) - T(i,j)) - (dt/hx^2)*(T_new(i-1,j) - 2*T_new(i,j) + T_new(i+1,j)) - (dt/hy^2)*(T_new(i,j-1) - 2*T_new(i,j) + T_new(i,j+1)))^2;
                    end
                end
                norm = sqrt(norm/(Nx*Ny));
            end
        end

        function tab = get_stability(T)
            stability = strings(size(T));
            for i=1:size(T,1)
                for j=1:size(T,2)
                    % If method is unstable, it will either give massive
                    % overshoots or simply NaN, if overshoot is higher than
                    % realmax
                    if(max(T{i,j},[],'all') >= 1. || sum(sum(isnan(T{i,j}))) > 0)
                        stability(i,j) = "Unstable";
                    else
                        stability(i,j) = "Stable";
                    end
                end
            end
            rowNames = {'dt = 1/64','dt = 1/128','dt = 1/256','dt = 1/512','dt = 1/1024','dt = 1/2048','dt = 1/4096'};
            varNames = {'Nx,Ny = 3','Nx,Ny = 7','Nx,Ny = 15','Nx,Ny = 31'};
            tab = array2table(categorical(stability)',"VariableNames",varNames,"RowNames",rowNames);
            
        end
  
    end
end