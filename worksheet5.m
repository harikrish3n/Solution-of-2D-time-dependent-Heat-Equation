%% Time-Dependant Partial Differential Equations (Group 4)
% The folder contains three files : worksheet5.m, Numerical_Methods.m and Utilities.m
% worksheet4.m is the skeleton that contains no function definitions.
% Numerical_Methods.m contains the functions for explicit and implicit
% Euler methods.
% Utilities.m contains additional functions for plotting, saving images, etc.

clear; close all;

%% Explicit Euler Method

%Setting up parameters
Nx=[3,7,15,31]; Ny=[3,7,15,31];
dt=[1/64,1/128,1/256,1/512,1/1024,1/2048,1/4096];
img=1; %To keep track of number of figures

% Creating cell array to store the solution for every time-step and grid
% spacing. Cell arrays are used here as the size of the temperature vector
% is different for every grid spacing
T_expl=cell(size(Nx,2),size(dt,2));

% Initializing temperature to 1 inside the domain and zero at the
% boundaries. Here, the size of the temperature matrix is (Nx+2,Ny+2) as we
% use a ghost layer approach to treat the boundaries. This makes the
% looping easier without having to do condition checks.
for i=1:size(T_expl,1)
    for j=1:size(T_expl,2)
    T_expl{i,j}=zeros(Nx(i)+2,Ny(i)+2);
    T_expl{i,j}(2:Nx(i)+1,2:Ny(i)+1) = 1;
    end
end

% Single time loop for all time-steps (Does include if-else statements)
% Can also do it in other ways, but there is no major difference in
% computational runtime and this is more convenient for generating plots
for time=dt(end):dt(end):4/8
    for i = 1:size(Nx,2)
        T_expl{i,7} = Numerical_Methods.expl_euler(Nx(i),Ny(i),dt(7),T_expl{i,7});
        if(mod(time,dt(1))==0)
            T_expl{i,1}=Numerical_Methods.expl_euler(Nx(i),Ny(i),dt(1),T_expl{i,1});
            T_expl{i,2}=Numerical_Methods.expl_euler(Nx(i),Ny(i),dt(2),T_expl{i,2});
            T_expl{i,3}=Numerical_Methods.expl_euler(Nx(i),Ny(i),dt(3),T_expl{i,3});
            T_expl{i,4}=Numerical_Methods.expl_euler(Nx(i),Ny(i),dt(4),T_expl{i,4});
            T_expl{i,5}=Numerical_Methods.expl_euler(Nx(i),Ny(i),dt(5),T_expl{i,5});
            T_expl{i,6}=Numerical_Methods.expl_euler(Nx(i),Ny(i),dt(6),T_expl{i,6});
        elseif(mod(time,dt(2))==0)
            T_expl{i,2}=Numerical_Methods.expl_euler(Nx(i),Ny(i),dt(2),T_expl{i,2});
            T_expl{i,3}=Numerical_Methods.expl_euler(Nx(i),Ny(i),dt(3),T_expl{i,3});
            T_expl{i,4}=Numerical_Methods.expl_euler(Nx(i),Ny(i),dt(4),T_expl{i,4});
            T_expl{i,5}=Numerical_Methods.expl_euler(Nx(i),Ny(i),dt(5),T_expl{i,5});
            T_expl{i,6}=Numerical_Methods.expl_euler(Nx(i),Ny(i),dt(6),T_expl{i,6});
        elseif(mod(time,dt(3))==0)
            T_expl{i,3}=Numerical_Methods.expl_euler(Nx(i),Ny(i),dt(3),T_expl{i,3});
            T_expl{i,4}=Numerical_Methods.expl_euler(Nx(i),Ny(i),dt(4),T_expl{i,4});
            T_expl{i,5}=Numerical_Methods.expl_euler(Nx(i),Ny(i),dt(5),T_expl{i,5});
            T_expl{i,6}=Numerical_Methods.expl_euler(Nx(i),Ny(i),dt(6),T_expl{i,6});
        elseif(mod(time,dt(4))==0)
            T_expl{i,4}=Numerical_Methods.expl_euler(Nx(i),Ny(i),dt(4),T_expl{i,4});
            T_expl{i,5}=Numerical_Methods.expl_euler(Nx(i),Ny(i),dt(5),T_expl{i,5});
            T_expl{i,6}=Numerical_Methods.expl_euler(Nx(i),Ny(i),dt(6),T_expl{i,6});
        elseif(mod(time,dt(5))==0)
            T_expl{i,5}=Numerical_Methods.expl_euler(Nx(i),Ny(i),dt(5),T_expl{i,5});
            T_expl{i,6}=Numerical_Methods.expl_euler(Nx(i),Ny(i),dt(6),T_expl{i,6});
        elseif(mod(time,dt(6))==0)
            T_expl{i,6}=Numerical_Methods.expl_euler(Nx(i),Ny(i),dt(6),T_expl{i,6});
        end
    end

    if(time== 1/8 || time == 2/8 || time== 3/8 || time== 4/8)
        figure(img);
        Utilities.surface_plot_expl(T_expl,Nx,Ny,dt,time); 
        Utilities.generate_plot(T_expl,Nx,Ny,dt,time);
        img=img+1;
    end
end

% Checking stability of explicit Euler method for every case
% MATLAB doesn't display more than 5 columns so the order of rows and
% columns is reversed here. Row -> dt, Column -> Nx,Ny
disp(table(Numerical_Methods.get_stability(T_expl),'VariableNames',{'Stability of explicit Euler method'}));

%% Implicit Euler Method

% Setting up parameters and cell array
dt_GS = 1/64;
T_impl=cell(size(Nx,2),1);

%Initialization
for i=1:size(T_impl,1)
    T_impl{i}=zeros(Nx(i)+2,Ny(i)+2);
    T_impl{i}(2:Nx(i)+1,2:Ny(i)+1) = 1;
end

%Time-loop
for time=dt_GS:dt_GS:4/8
    for i = 1:size(Nx,2)
        T_impl{i} = Numerical_Methods.impl_euler(Nx(i),Ny(i),dt_GS,T_impl{i});
    end
    if(time== 1/8 || time == 2/8 || time== 3/8 || time== 4/8)
        figure(img);
        Utilities.surface_plot_impl(T_impl,Nx,Ny,dt_GS,time);        
        img=img+1;
    end
end