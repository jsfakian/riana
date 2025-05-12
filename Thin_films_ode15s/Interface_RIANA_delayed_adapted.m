% Usage from MATLAB:
%   result = calc(Ep1, wavelength1, tp1, t_delay1, t_max1, L1, ...
%                   material, material_substrate, n1, k1, n2, k2);
%
% Inputs:
%   Ep1                - pulse fluence (J/cm^2)
%   wavelength1        - wavelength (nm)
%   tp1                - pulse duration (fs)
%   t_delay1           - pulse separation / delay (fs)
%   t_max1             - maximum time (ps)
%   L1                 - film thickness (nm)
%   material           - one of {'Cr','Ti','Ni','Au','Ag','Al','Cu','W','Pt','Steel','Mo'}
%   material_substrate - one of {'SiO2','Si','Soda lime glass','Air'}
%   n1, k1             - optical constants for the film
%   n2, k2             - optical constants for the substrate
%
% Output:
%   result             - whatever your Code_for_* routines return

function result = calc( ...
    Ep1, wavelength1, tp1, t_delay1, t_max1, L1, ...
    material, material_substrate, n1, k1, n2, k2)

    % Suppress warnings and plots by default
    warning('off');
    global stri
    stri = 'off';

    % Assign globals for downstream code
    global n_1 k_1 n_2 k_2
    n_1 = n1;
    k_1 = k1;
    n_2 = n2;
    k_2 = k2;

    %% SUBSTRATE
    if strcmp(material_substrate, 'Si') == 1
        if wavelength1 / 1000 == 0.800
            n_2 = 3.6750; % Green
            k_2 = 0.0054113;

        elseif wavelength1 / 1000 == 1.026
            n_2 = 3.5632; % Green
            k_2 = 0.00027806;

        elseif wavelength1 / 1000 == 0.515
            n_2 = 4.2170; % Green
            k_2 = 0.037;

        elseif wavelength1 / 1000 == 0.248
            n_2 = 1.57; % Aspnes
            k_2 = 3.5650;
        end

    elseif strcmp(material_substrate, 'Soda') == 1
        n_2 = 1.5130 - 0.003169 * (wavelength1 / 1000)^2 + 0.003962 / (wavelength1 / 1000)^2;
        k_2 = 0;

    elseif strcmp(material_substrate, 'SiO2') == 1
        wavelength = wavelength1 * 1e-3; % in mum
        term1 = 0.6961663 * wavelength^2 / (wavelength^2 - 0.0684043^2);
        term2 = 0.4079426 * wavelength^2 / (wavelength^2 - 0.1162414^2);
        term3 = 0.8974794 * wavelength^2 / (wavelength^2 - 9.896161^2);
        n_2 = sqrt(1 + term1 + term2 + term3);
        k_2 = 0;

    elseif strcmp(material_substrate, 'Air') == 1
        n_2 = 1;
        k_2 = 0;
    end

    %% UPPER MATERIAL (METAL)
    if strcmp(material, 'Mo') == 1 && wavelength1 == 1026
        n_1 = 2.4357;
        k_1 = 4.1672;
    end

    %% RUN THE CODES
    % Determine which backend to use
    isBigDelay = (t_delay1 > 5000);
    isTwoMetalRoutine = ...
        any(strcmp(material, {'Ni','Steel','Ti','Pt'})) || ...
        ( any(strcmp(material, {'Au','Cu','Cr','Ag','W','Al','Mo'})) && isBigDelay );

    % Dispatch
    if any(strcmp(material, {'Au','Cu','Cr','Ag','W','Al','Mo'})) && ~isBigDelay
        % Use the single-metal backend
        result = Code_for_all_metals_but_two( ...
            Ep1, wavelength1, tp1, t_delay1, t_max1, ...
            material, material_substrate, L1);
    elseif isTwoMetalRoutine
        % Use the two-metal / delayed backend
        result = Code_for_two( ...
            Ep1, wavelength1, tp1, t_delay1, t_max1, ...
            material, material_substrate, L1);
    else
        % No matching routine
        warning('No simulation routine for material=%s, delay=%g fs', material, t_delay1);
        result = [];
    end
end
