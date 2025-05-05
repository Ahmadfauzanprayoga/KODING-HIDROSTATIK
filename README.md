%KODING BUATAN AHMAD FAUZAN PRAYOGA DILARANG DIBAGIKAN TANPA SEPENGETAHUAN PEMILIK UNTUK MATHLAB
clc;
clear all;
close all;

% -----------------------------
% INPUT DATA FROM EXCEL
% SKRIP JADIKAN SATU FOLDER DENGAN EXCEL
nama_file = 'HALFBREADTH.xlsx';%MASUKAN NAMA SHEET EXCEL
sheet = 'Sheet1'; %MASUKAN NAMA SHEET EXCEL
range_data = 'C2:P24';  % RANGE EXCEL HALFBREADTH
[~, ~, raw_data] = xlsread(nama_file, sheet, range_data);
half_breadth = zeros(size(raw_data)) * NaN;

for i = 1:size(raw_data, 1)
    for j = 1:size(raw_data, 2)
        if isnumeric(raw_data{i,j})
            half_breadth(i,j) = raw_data{i,j};
        elseif strcmp(raw_data{i,j}, '-')
            half_breadth(i,j) = NaN;
        elseif ischar(raw_data{i,j}) && ~isempty(str2double(raw_data{i,j}))
            half_breadth(i,j) = str2double(raw_data{i,j});
        else
            half_breadth(i,j) = NaN;
        end
    end
end
half_breadth(isnan(half_breadth)) = 0;

% -----------------------------
% PARAMETER KAPAL
% -----------------------------
delta_L = 1;           % MASUKKAN Station spacing (m)
delta_WL = 1;          % MASUKKAN Waterline spacing (m)
rho = 1025;            % MASUKKAN Seawater density (kg/m³)

[n_stations, n_waterlines] = size(half_breadth);
stations = (0:n_stations-1) * delta_L;
midship_pos = mean(stations);
Lpp = max(stations);
drafts = (0:n_waterlines-1) * delta_WL;

% Station integration coefficients
if mod(n_stations, 2) == 1 % Simpson's rule (odd stations)
    coeff_stations = ones(1, n_stations);
    coeff_stations(2:2:end-1) = 4;
    coeff_stations(3:2:end-2) = 2;
    delta_station = delta_L / 3;
else % Trapezoidal rule (even stations)
    coeff_stations = 2 * ones(1, n_stations);
    coeff_stations([1 end]) = 1;
    delta_station = delta_L / 2;
end

% -----------------------------
% VARIABEL HASIL YANG AKAN KELUAR
% -----------------------------
all_drafts = [];
all_LCB = [];
all_LCF = [];
all_Cb = [];
all_Cm = [];
all_Cp = [];
all_Cw = [];
all_Volume = [];
all_Displacement = [];
all_KB = [];
all_BM = [];
all_KM = [];
all_WSA = [];
all_TPC = [];
all_MCT = [];
all_Am = [];

% -----------------------------
% PERHITUNGAN PER WATERLINE
% -----------------------------
for wl = 1:n_waterlines
    T = drafts(wl);
    if T == 0, continue; end % Skip zero draft
    
    % 1. Sectional Area, Volume, Displacement, KB
    A_station = zeros(1, n_stations);
    momen_KB = zeros(1, n_stations);
    for i = 1:n_stations
        y = half_breadth(i, 1:wl);
        n = length(y);
        if n < 2
            A_station(i) = 0;
            momen_KB(i) = 0;
            continue;
        end
        
        % Waterline integration coefficients
        if mod(n, 2) == 1 % Simpson's
            coeff = ones(1, n);
            coeff(2:2:end-1) = 4;
            coeff(3:2:end-2) = 2;
            delta_waterline = delta_WL / 3;
        else % Trapezoidal
            coeff = 2 * ones(1, n);
            coeff([1 end]) = 1;
            delta_waterline = delta_WL / 2;
        end
        
        A_station(i) = delta_waterline * sum(coeff .* y) * 2; % Both sides
        
        % Moment for KB
        z_vals = drafts(1:wl);
        sum_z = sum(coeff .* z_vals);
        z_centroid = sum_z / sum(coeff);
        momen_KB(i) = A_station(i) * z_centroid;
    end
    
    volume = delta_station * sum(coeff_stations .* A_station);
    displacement = volume * rho;
    total_momen_KB = delta_station * sum(coeff_stations .* momen_KB);
    KB = total_momen_KB / volume;
    
    % 2. Wetted Surface Area (WSA)
    G_station = zeros(1, n_stations);
    for i = 1:n_stations
        y = half_breadth(i, 1:wl);
        if length(y) < 2
            G_station(i) = 0;
            continue;
        end
        g = sum(sqrt(delta_WL^2 + diff(y).^2));
        G_station(i) = 2 * g; % Both sides
    end
    WSA = delta_station * sum(coeff_stations .* G_station);
    
    % 3. LCB & LCF from Midship
    y_waterline = half_breadth(:, wl);
    luas_garis_air = 2 * delta_station * sum(coeff_stations .* y_waterline');
    
    momen_LCB = A_station .* stations;
    total_momen_LCB = delta_station * sum(coeff_stations .* momen_LCB);
    LCB = total_momen_LCB / volume;
    
    momen_LCF = y_waterline .* stations';
    total_momen_LCF = 2 * delta_station * sum(coeff_stations .* momen_LCF');
    LCF = total_momen_LCF / luas_garis_air;
    
    % 4. KM & BM
    I_T = delta_station * sum(coeff_stations .* ((2/3) * y_waterline.^3)');
    BM = I_T / volume;
    KM = KB + BM;
    
    % 5. TPC
    TPC = (luas_garis_air * rho) / 100;
    
    % 6. MCT
    x_midship = stations - midship_pos;
    I_L = delta_station * sum(coeff_stations .* (2 * y_waterline .* x_midship'.^2)');
    BML = I_L / volume;
    MCT = (displacement * BML) / Lpp;
    
    % 7. Midship Section & Coefficients
    [~, idx_mid] = min(abs(stations - midship_pos));
    y_mid = half_breadth(idx_mid, 1:wl);
    n_mid = length(y_mid);
    if n_mid >= 2
        if mod(n_mid, 2) == 1 % Simpson's
            coeff_mid = [1, repmat([4, 2], 1, (n_mid-3)/2), 4, 1];
            delta_mid = delta_WL / 3;
        else % Trapezoidal
            coeff_mid = [1, 2*ones(1, n_mid-2), 1];
            delta_mid = delta_WL / 2;
        end
        Am = 2 * delta_mid * sum(coeff_mid .* y_mid);
    else
        Am = 0;
    end
    
    B = 2 * max(y_waterline);
    Cb = volume / (Lpp * B * T);
    Cm = Am / (B * T);
    Cp = volume / (Am * Lpp);
    Cw = luas_garis_air / (Lpp * B);
    
    % Store results
    all_drafts(end+1) = T;
    all_LCB(end+1) = LCB - midship_pos;
    all_LCF(end+1) = LCF - midship_pos;
    all_Cb(end+1) = Cb;
    all_Cm(end+1) = Cm;
    all_Cp(end+1) = Cp;
    all_Cw(end+1) = Cw;
    all_Volume(end+1) = volume;
    all_Displacement(end+1) = displacement / 1000;
    all_KB(end+1) = KB;
    all_BM(end+1) = BM;
    all_KM(end+1) = KM;
    all_WSA(end+1) = WSA;
    all_TPC(end+1) = TPC;
    all_MCT(end+1) = MCT;
    all_Am(end+1) = Am;
    
    % Display
    fprintf('=== Draft: %.2f m ===\n', T);
    fprintf('Volume      : %.2f m³\n', volume);
    fprintf('Displacement: %.2f ton\n', displacement/1000);
    fprintf('KB          : %.2f m\n', KB);
    fprintf('BM          : %.2f m\n', BM);
    fprintf('KM          : %.2f m\n', KM);
    fprintf('LCB         : %.2f m\n', all_LCB(end));
    fprintf('LCF         : %.2f m\n', all_LCF(end));
    fprintf('WSA         : %.2f m²\n', WSA);
    fprintf('TPC         : %.2f ton/cm\n', TPC);
    fprintf('MCT         : %.2f ton·m/cm\n', MCT);
    fprintf('Amidship    : %.2f m²\n', Am);
    fprintf('Cb: %.3f, Cm: %.3f, Cp: %.3f, Cw: %.3f\n', Cb, Cm, Cp, Cw);
    fprintf('-----------------------------\n');
end

% -----------------------------
% PLOTTING
% -----------------------------
% Grafik LCB vs LCF
figure('Name','LCB & LCF vs Draft','NumberTitle','off');
plot(all_LCB, all_drafts, 'b-o', all_LCF, all_drafts, 'r-s');
ylabel('Draft (m)');
xlabel('Posisi dari Midship (m)');
legend('LCB', 'LCF', 'Location','best');
title('LCB dan LCF vs Draft');
grid on;
set(gca, 'YDir','normal');

% Grafik Koefisien
figure('Name','Koefisien vs Draft','NumberTitle','off');
hold on;
plot(all_Cb, all_drafts, 'b-o', 'LineWidth',1.5, 'DisplayName','Cb');
plot(all_Cm, all_drafts, 'r--s', 'LineWidth',1.5, 'DisplayName','Cm');
plot(all_Cp, all_drafts, 'g-.^', 'LineWidth',1.5, 'DisplayName','Cp');
plot(all_Cw, all_drafts, 'm:d', 'LineWidth',1.5, 'DisplayName','Cw');
hold off;
ylabel('Draft (m)');
xlabel('Nilai Koefisien');
legend('show','Location','best');
title('Koefisien Bentuk Kapal vs Draft');
grid on;
set(gca, 'YDir','normal');

% Grafik Individual
parameters = {
    {'Volume', all_Volume, 'Volume (m³)'},...
    {'Displacement', all_Displacement, 'Displacement (ton)'},...
    {'KB', all_KB, 'KB (m)'},...
    {'BM', all_BM, 'BM (m)'},...
    {'KM', all_KM, 'KM (m)'},...
    {'WSA', all_WSA, 'WSA (m²)'},...
    {'TPC', all_TPC, 'TPC (ton/cm)'},...
    {'MCT', all_MCT, 'MCT (ton·m/cm)'},...
    {'Amidship Area', all_Am, 'Amidship Area (m²)'}
};

for i = 1:length(parameters)
    figure('Name',parameters{i}{1},'NumberTitle','off');
    plot(parameters{i}{2}, all_drafts, 'k-o', 'LineWidth',1.5);
    ylabel('Draft (m)');
    xlabel(parameters{i}{3});
    title([parameters{i}{1} ' vs Draft']);
    grid on;
    set(gca, 'YDir','normal');
end
