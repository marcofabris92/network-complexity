close all
clearvars
clc

%% Loading data and computing index Phi_G in the deterministic setting

% CElegans
[CEw,CEb,clustersCE,KCE] = CElegansData();
Phi_G_CEw_DET = netcompl_(CEw,KCE,clustersCE);
Phi_G_CEb_DET = netcompl_(CEb,KCE,clustersCE);
Phi_G_CEw_DET
Phi_G_CEb_DET
Phi_G_CEw_DET_ = Phi_G_CEw_DET/numnodes(CEw);
Phi_G_CEb_DET_ = Phi_G_CEb_DET/numnodes(CEb);
Phi_G_CEw_DET_
Phi_G_CEb_DET_

% Power Grid
[PG,clustersPG,KPG] = PowerGridData();
Phi_G_PG_DET = netcompl_(PG,KPG,clustersPG);
Phi_G_PG_DET
Phi_G_PG_DET_ = Phi_G_PG_DET/numnodes(PG);
Phi_G_PG_DET_

% Opinion Dynamics
[OD,clustersOD,KOD] = OpinionDynamicsData();
Phi_G_OD_DET = netcompl_(OD,KOD,clustersOD);
Phi_G_OD_DET
Phi_G_OD_DET_ = Phi_G_OD_DET/numnodes(OD);
Phi_G_OD_DET_

%% Monte Carlo runner
NMC = 1e2;      % number of runs
type = 1;       % rewiring type {0,1,2,3}
p = 1;          % rewiring prob. (only for types 0 and 3; default: p=1)

% CElegans
par.n = numnodes(CEb);
par.m = numedges(CEb);
par.weights = reshape(adjacency(CEw),par.n^2,1);
par.weights = par.weights(par.weights~=0);
par.p = p;
par.Q = 1:(par.n^2 - par.n);
par.clusters = clustersCE;
par.K = KCE;
Phi_Gs_CEw = zeros(NMC,1);      % compl. indices for CEw
Phi_Gs_CEb = zeros(NMC,1);      % compl. indices for CEb
parfor j = 1:NMC
    Phi_Gs_CEw(j) = netcompl_(rewiring(CEw,par,type),KCE,clustersCE);
    Phi_Gs_CEb(j) = netcompl_(rewiring(CEb,par,type),KCE,clustersCE);
end

% Power Grid
par.n = numnodes(PG);
par.m = numedges(PG);
par.weights = reshape(adjacency(PG),par.n^2,1);
par.weights = par.weights(par.weights~=0);
par.p = p;
par.Q = 1:(par.n^2 - par.n);
par.clusters = clustersPG;
par.K = KPG;
Phi_Gs_PG = zeros(NMC,1);      % compl. indices for PG
parfor j = 1:NMC
    Phi_Gs_PG(j) = netcompl_(rewiring(PG,par,type),KPG,clustersPG);
end


% Opinion Dynamics
par.n = numnodes(OD);
par.m = numedges(OD);
par.weights = reshape(adjacency(OD),par.n^2,1);
par.weights = par.weights(par.weights~=0);
par.p = p;
par.Q = 1:(par.n^2 - par.n);
par.clusters = clustersOD;
par.K = KOD;
Phi_Gs_OD = zeros(NMC,1);      % compl. indices for OD
parfor j = 1:NMC
    Phi_Gs_OD(j) = netcompl_(rewiring(OD,par,type),KOD,clustersOD);
end

Phi_Gs_CEw
Phi_Gs_CEb
Phi_Gs_PG
Phi_Gs_OD


%% Computing mean and std. dev.

% mean
mu_CEw = mean(Phi_Gs_CEw)/numnodes(CEw);
mu_CEb = mean(Phi_Gs_CEb)/numnodes(CEb);
mu_PG = mean(Phi_Gs_PG)/numnodes(PG);
mu_OD = mean(Phi_Gs_OD)/numnodes(OD);
mu_CEw
mu_CEb
mu_PG
mu_OD

% std. dev.
std_CEw = std(Phi_Gs_CEw)/numnodes(CEw);
std_CEb = std(Phi_Gs_CEb)/numnodes(CEb);
std_PG = std(Phi_Gs_PG)/numnodes(PG);
std_OD = std(Phi_Gs_OD)/numnodes(OD);
std_CEw
std_CEb
std_PG
std_OD




%% END OF MAIN %%


%%%%%%%%%%%%%%%%%%
%%%    DATA    %%%
%%%%%%%%%%%%%%%%%%



%% CElegans data
function [CEw,CEb,clustersCE,KCE] = CElegansData()

% Importing the deterministic graph under analysis
[A,Neuron_ordered,class] = datareader('chem','weighted');

% load the GABAergic synapse listing
GABAergic = GABA; 

% treat GABAergic synapses in chemical network as inhibitory
A(find(GABAergic),:) = -A(find(GABAergic),:);

% number of nodes
n = length(A);

% number of clusters imposed
KCE = 3;
if n < KCE 
    error(['The network must have at most ' num2str(n) ' clusters.'])
end

% determine the category of neuron, i.e. the clusters
sen = [];
int = [];
mot = [];

parfor i = 1:n
    if (findstr(char(class(i)),'S') > 1)
        sen = [sen i];
    elseif (findstr(char(class(i)),'M') > 1)
        mot = [mot i];
    elseif (findstr(char(class(i)),'I') > 1)
        int = [int i];
    end
end

categories = cell(KCE,1);
categories{1} = sen;
categories{2} = int;
categories{3} = mot;

clustersCE = zeros(n,1);
for k = 1:KCE
    for i = 1:length(categories{k})
        clustersCE(categories{k}(i)) = k;
    end
end


% the three cluster sizes
% n
% sen_length = length(sen)
% int_length = length(int)
% mot_length = length(mot)

% reshape matrix so partitioned by category, anteroposterior order within
AA = [A(sen,sen) A(sen,int) A(sen,mot);
      A(int,sen) A(int,int) A(int,mot);
      A(mot,sen) A(mot,int) A(mot,mot)];

for i = 1:n
    AA(i,i) = 0;
end

% CEw: CElegans weighted graph; CEb: CEw but weights are binarized
CEw = digraph(full(AA));
CEb = full(AA);
CEb(CEb~=0) = 1;
CEb = digraph(CEb);

end


%% Power Grid data
function [PG,clustersPG,KPG] = PowerGridData()

EdgeTable = table(1+[0 1; 1 2; 2 3; 3 4; 4 5; 5 6; 6 7; 7 0; 7 40; 0 8; 8 5; 8 3; 2 9; 9 10; 10 11; 11 12; 12 3; 11 9; 8 13; 9 14; 14 15; 15 16; 16 17; 17 18; 18 15; 14 19; 19 20; 20 21; 21 15; 21 22; 22 23; 23 24; 24 21; 23 25; 25 26; 26 27; 27 28; 28 29; 29 26; 29 30; 30 31; 31 25; 28 32; 32 30; 32 33; 33 34; 33 35; 35 32; 33 36; 36 32; 32 37; 37 38; 38 36; 39 40; 40 31; 31 39; 40 41; 41 37; 41 42; 42 43; 43 41; 33 44; 44 45; 45 46; 46 44; 44 47; 47 48; 47 223; 48 36; 33 49; 49 44; 49 50; 50 51; 51 38; 51 52; 52 53; 53 54; 55 53; 53 56; 56 57; 57 58; 58 59; 58 60; 60 61; 61 62; 61 63; 61 64; 53 65; 65 66; 66 50; 48 67; 67 68; 68 69; 69 70; 70 57; 57 69; 69 71; 71 72; 70 71; 71 73; 73 74; 74 75; 75 76; 76 77; 77 78; 78 79; 79 80; 80 81; 81 82; 82 83; 83 84; 84 71; 84 112; 43 85; 85 57; 85 90; 86 87; 87 88; 88 89; 88 90; 90 91; 91 92; 92 93; 93 38; 92 65; 65 94; 94 95; 95 96; 96 97; 97 98; 97 99; 96 100; 100 101; 101 102; 96 103; 103 104; 104 105; 105 106; 106 107; 57 37; 37 96; 53 37; 67 108; 108 48; 108 109; 109 110; 110 111; 111 68; 68 96; 111 112; 112 113; 113 110; 113 114; 113 115; 115 67; 57 60; 60 116; 116 117; 117 118; 118 119; 119 120; 119 76; 117 121; 121 122; 122 123; 122 124; 124 125; 125 126; 126 127; 125 128; 128 129; 129 130; 130 131; 131 117; 122 132; 132 133; 133 134; 134 135; 135 136; 135 137; 135 138; 138 139; 139 140; 140 141; 141 142; 142 143; 142 144; 144 140; 144 145; 145 146; 146 147; 147 148; 148 149; 149 150; 149 151; 148 152; 152 153; 153 154; 154 155; 155 156; 156 157; 157 158; 158 159; 159 160; 160 13; 13 161; 161 157; 161 162; 162 163; 163 13; 163 164; 164 165; 165 154; 165 147; 164 166; 166 167; 167 168; 168 42; 42 162; 42 169; 169 170; 170 152; 162 171; 171 153; 153 169; 168 166; 167 125; 125 172; 172 167; 167 173; 173 166; 166 174; 174 173; 173 175; 175 174; 174 176; 176 177; 177 175; 175 178; 178 177; 179 180; 180 181; 181 182; 182 173; 181 175; 111 183; 183 184; 184 185; 184 109; 184 186; 186 187; 187 183; 187 188; 187 189; 187 190; 190 191; 191 192; 192 193; 191 194; 194 195; 195 196; 196 197; 195 198; 198 197; 197 199; 199 200; 200 201; 201 202; 200 203; 203 202; 202 204; 204 81; 190 183; 183 205; 186 183; 186 206; 206 207; 207 208; 208 209; 208 210; 210 211; 211 206; 210 212; 212 186; 210 213; 213 212; 213 214; 214 186; 214 215; 215 216; 216 217; 217 218; 218 219; 219 220; 220 221; 221 222; 222 223; 223 224; 224 225; 225 219; 225 226; 226 227; 227 220; 227 228; 228 217; 225 229; 229 224; 229 226; 229 230; 229 231; 231 232; 232 233; 233 212; 233 229; 232 234; 234 233; 234 186; 212 226; 226 235; 235 213], 'VariableNames', {'EndNodes'});

PG = digraph(EdgeTable);

AA = adjacency(PG);

for i = 1:numnodes(PG)
    AA(i,i) = 0;
end

PG = digraph(AA);

clustersPG = 1+[0 0 1 1 1 0 1 0 0 0 1 0 1 1 1 1 1 1 0 0 0 0 0 1 0 1 0 0 0 0 1 0 1 1 1 1 0 0 0 0 0 0 0 0 1 1 0 1 0 0 0 1 1 1 1 0 0 0 0 0 0 1 1 1 1 0 0 1 0 1 0 0 0 1 1 0 1 0 0 0 0 1 1 1 0 1 0 1 0 1 0 1 1 0 0 1 0 1 0 0 0 0 0 1 1 1 1 0 0 1 1 1 0 0 0 1 1 0 0 1 0 1 0 1 1 1 0 1 0 1 1 0 0 0 1 0 0 1 1 1 1 0 0 1 1 0 0 1 0 1 1 0 0 1 0 0 0 1 1 0 1 1 0 1 1 1 1 1 1 0 0 0 0 1 0 0 1 1 0 0 0 0 1 0 1 1 1 1 1 1 1 0 1 0 0 1 0 0 1 1 1 0 0 1 1 1 0 1 1 0 0 0 1 1 0 1 0 1 1 0 1 1 0 1 0 1 0 1 0 0 1 1 1 0 1 1];

KPG = 2;

end


%% Opinion Dynamics data
function [OD,clustersOD,KOD] = OpinionDynamicsData()

% edge list info
fileName = 'polblogs_edges.gml';
inputfile = fopen(fileName);
l = 0;
k = 1;
EdgeTable = zeros([],2) ;
while true
    % Get a line from the input file
    tline = fgetl(inputfile);
    % Quit if end of file
    if ~ischar(tline)
        break
    end
    nums = regexp(tline,'\d+','match');
    if ~isempty(nums)
        if l == 1
            l = 0;
            EdgeTable(k,2) = str2double(nums{1});
            k = k + 1;
            continue
        end
        EdgeTable(k,1) = str2double(nums{1});
        l = 1;
    else
        l = 0;
        continue
    end
end
fclose(inputfile);

EdgeTable = table(EdgeTable, 'VariableNames', {'EndNodes'});

OD = digraph(EdgeTable);

AA = adjacency(OD);

for i = 1:numnodes(OD)
    AA(i,i) = 0;
end

OD = digraph(AA);


% nodes' clusters info
fileName = 'polblogs_nodes.gml';
inputfile = fopen(fileName);
clustersOD = zeros([],1) ;
id = [];
ID = 0;
value = [];
while true
    % Get a line from the input file
    tline = fgetl(inputfile);

    % Quit if end of file
    if ~ischar(tline)
        break
    end

    id = regexp(tline,'id ') + 3;
    if ~isempty(id)
        I = id;
        id = [];
        for i = I:length(tline)
            id = [id char(extract(tline,i))];
        end
        ID = str2double(id);
    end
    
    value = regexp(tline,'value ') + 6;
    if ~isempty(value)
        I = value;
        value = [];
        for i = I:length(tline)
            value = [value char(extract(tline,i))];
        end
        clustersOD(ID) = 1+str2double(value);
    end

end
fclose(inputfile);

KOD = 2;


end
