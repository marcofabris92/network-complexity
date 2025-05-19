clc 
clear all
close all

n = 279;

[A,Neuron_ordered,class] = datareader('chem','weighted');

%load the GABAergic synapse listing
GABAergic = GABA; 

%treat GABAergic synapses in chemical network as inhibitory
A(find(GABAergic),:) = -A(find(GABAergic),:);

%number of nodes
n = length(A);

%determine the category of neuron
sen = [];
int = [];
mot = [];

for ii = 1:n
    if (findstr(char(class(ii)),'S') > 1)
        sen = [sen ii];
    elseif (findstr(char(class(ii)),'M') > 1)
        mot = [mot ii];
    elseif (findstr(char(class(ii)),'I') > 1)
        int = [int ii];
    end
end

length(sen)
length(int)
length(mot)

%reshape matrix so partitioned by category, anteroposterior order within
AA = [A(sen,sen) A(sen,int) A(sen,mot);
      A(int,sen) A(int,int) A(int,mot);
      A(mot,sen) A(mot,int) A(mot,mot)];

A = full(AA);
A_bin = full(AA);
A_bin(A_bin~=0)=1;


figure;
imagesc(A_bin);
colormap(gray);
axis equal tight;
title('Binarized Adjacency Matrix (Grouped by Neuron Type)');
xlabel('Neuron Index');
ylabel('Neuron Index');

% Add partition lines
hold on;
% Indices where groups end
idx_sen = length(sen);
idx_int = idx_sen + length(int);
idx_mot = idx_int + length(mot);

% Draw horizontal and vertical lines at boundaries
xline(idx_sen + 0.5, 'r', 'LineWidth', 1.2); % between sensory and interneuron
xline(idx_int + 0.5, 'r', 'LineWidth', 1.2); % between interneuron and motor
yline(idx_sen + 0.5, 'r', 'LineWidth', 1.2);
yline(idx_int + 0.5, 'r', 'LineWidth', 1.2);

% Optional: add text labels
text(idx_sen/2, -10, 'Sensory', 'HorizontalAlignment', 'center', 'Color', 'b');
text((idx_sen + idx_int)/2, -10, 'Interneuron', 'HorizontalAlignment', 'center', 'Color', 'b');
text((idx_int + idx_mot)/2, -10, 'Motor', 'HorizontalAlignment', 'center', 'Color', 'b');
