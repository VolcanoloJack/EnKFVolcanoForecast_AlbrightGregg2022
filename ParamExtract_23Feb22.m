function [ val, vallo, valhi, cpu ] = ParamExtract_23Feb22( filenames )
% Given a set of synthetic model run filenames, extracts final parameters,
% calculates misfits, and normalizes to synthetic values

%% Changelog

% 28 Jun 21 - Original version based on ParamExtract_07Dec20.m

%% Setup

% Percentiles to be extracted
lop = 5; hip = 95;

% Number of model runs
if isa(filenames,'char')
    filenum = 1;
else
    filenum = length(filenames);
end

% Setup for extraction
w = 0;
val = zeros(filenum+1,10);
vallo = zeros(filenum+1,10);
valhi = zeros(filenum+1,10);
cpu = zeros(filenum,1);

misfits_tot = [];

%% Load Model Runs and Extract Data
for j = 1:filenum
    if isa(filenames,'char')
        load(filenames);
    else
        load(filenames(j));
    end
    w = w+1;
    
    cpu(w) = cpusteps(end);
    
    misfits = zeros(13,N);
    
    % RMSE
    misfits(1,:) = params(end,:,end);
    
    % Tensile Stress
    k = length(steps_s);
    TSs = max(TFs{k,1}(1,:));
    
    for i = 1:N
        misfits(2,i) = (max(TF{i,end}(1,:)) - TSs)/abs(TSs);
    end
    
    % Unique Parameters [dV, aspect ratio, distance]
    if drv == 7
        misfits(3,:) = (params(1,:,end)-synths(1,end))/abs(synths(1,end));
        U = [5,6,12];
    else
        misfits(3,:) = (params(7,:,end)-synths(7,end))/abs(synths(7,end));
        U = [4,5,6];
    end
    
    misfits(4,:) = log(params(8,:,end))-log(synths(8,end));
    
    errdist = zeros(1,N);
    for i  = U
        errdist = errdist + (params(i,:,end)-synths(i,end)).^2;
    end
    misfits(5,:) = sqrt(errdist)/synths(U(3),end);
    
    misfits(6,:) = mean(abs(misfits(3:5,:)));
    
    % Nonunique Parameters [r1, r2, dP]
    if drv == 7
        U = 9:10;
        misfits(9,:) = (params(11,:,end)-synths(11,end))/50e6;
    else
        U = 1:3;
        misfits(9,:) = ((params(3,:,end)*pscale(3))-synths(3,end))/50e6;
    end
    
    for i = 1:2
        k = U(i);
        misfits(6+i,:) = (params(k,:,end)-synths(k,end))/abs(synths(k,end));
    end
    
    misfits(10,:) = mean(abs(misfits(7:9,:)));
    
    % Process terms
    for i = 1:10
        val(w,i) = mean(misfits(i,:));
        vallo(w,i) = abs(mean(misfits(i,:)) - prctile(misfits(i,:),lop));
        valhi(w,i) = abs(mean(misfits(i,:)) - prctile(misfits(i,:),hip));
    end
    
    misfits_tot = [misfits_tot misfits];
    
    clearvars -except filenames val vallo valhi w filenum j lop hip...
        misfits_tot cpu
end

for i = 1:10
    val(end,i) = mean(misfits_tot(i,:));
    vallo(end,i) = abs(mean(misfits_tot(i,:)) - prctile(misfits_tot(i,:),lop));
    valhi(end,i) = abs(mean(misfits_tot(i,:)) - prctile(misfits_tot(i,:),hip));
end

end

