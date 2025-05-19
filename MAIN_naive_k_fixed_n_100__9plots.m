%%MAIN_naive_k
% fixed k, clustering: naive, N=100, tests=100

% clearing options
close all
clear all
%warning('off','all')
clc

% SAVE = 1 OR LOAD = 0 GRAPHS 
W = 1;

vrb = 0;
fprintf('Initializing\n\n')

% graphic options
pp.lw = 2;
pp.ftsz = 20;
pp.inter = 'latex';
pp.marksz = 5;

% simulation setup
pp.N = 100; %10:10:200;
pp.ks = [2 4]; % <---------------------------------------- k
tests = 1*1e2;
type = 'naive';
mymod = 0;
keys_ = [];

% clustering types available
types = {
    'mymodularity',...      % -2    BEGIN chosen modularity (specify mymod)
    'minmodularity',...     % -1          MIN modularity (mymod = -1/2)
    'maxmodularity',...     % 0     END   MAX modularity (mymod = 1)
    ...
    'indegree',...          % 1     BEGIN centrality
    'outdegree',...
    'incloseness',...
    'outcloseness',...
    'betweenness',...
    'pagerank',...
    'hubs',...
    'authorities',...       % 8     END centrality
    ...
    'neutral',...           % 9     BEGIN spectral
    'neutralsquared',...
    'symmetrized',...
    'symmetrizedsquared',...
    'inward',...
    'outward',...
    'inoutward',...
    'outinward',...
    'strong',...            % 17    END spectral    
    ...
    'random'...             % 18    BEGIN randomic
    'inrandom',...
    'outrandom',...
    'neutralrandom',...
    'symmetrizedrandom',...
    'purelyrandom',...
    'naive'};               % 24    END randomic
% IMPORTANT: You must fill key_trans with all the END-keys > 0!
key_trans = [8 17 24];
% IMPORTANT: You must specify "shift" as the number of nonpositive keys.
shift = 1+2;
shifted_dim = length(types)-shift;
dic = containers.Map(num2cell((1-shift:shifted_dim)),types);
dic_info = [mymod shift shifted_dim key_trans];
%for k = pp.ks %key = keys_
    
    %type = types{key};
    nodes = length(pp.N);
    clusters = length(pp.ks);
    Phi_Gs_ER1 = zeros(tests,clusters);
    Phi_Gs_ER2 = zeros(tests,clusters);
    Phi_Gs_ER3 = zeros(tests,clusters);
    Phi_Gs_SW1 = zeros(tests,clusters);
    Phi_Gs_SW2 = zeros(tests,clusters);
    Phi_Gs_SW3 = zeros(tests,clusters);
    Phi_Gs_BA1 = zeros(tests,clusters);
    Phi_Gs_BA2 = zeros(tests,clusters);
    Phi_Gs_BA3 = zeros(tests,clusters);
    avg_time_per_test = 0;
    for t = 1:tests
        tic
        
        % Gneration + Clustering phase
        for n_ = 1:length(pp.N)
            
            
            n = pp.N(n_);
            %% Generation of graphs
            % Random graph (Erdos-Renyi)
            %fprintf('ER generation\n')
            %pER = rand;
            ERthr = log(n)/n;
            pp.pER1 = 0.25*ERthr;
            pp.pER2 = ERthr;
            pp.pER3 = 1.25*ERthr;
            ER1 = ErdosRenyi(n,pp.pER1); 
            ER2 = ErdosRenyi(n,pp.pER2);
            ER3 = ErdosRenyi(n,pp.pER3);

            % Small-World (Watts-Strogatz)
            %fprintf('SW generation\n')
            %etaSW = rand;
            %K = floor((1+log(n)/2)*etaSW+(floor(n/2)-1)*(1-etaSW));
            K = 2;
            %pSW = rand;
            pSW = 0.5;
            %C_G_min = rand;
            C_G_min = 0;
            %connSW = round(rand);
            connSW = 0;
            SW1 = SmallWorld(n,K,pSW-0.25,C_G_min,connSW,vrb);
            SW2 = SmallWorld(n,K,pSW,C_G_min,connSW,vrb);
            SW3 = SmallWorld(n,K,pSW+0.25,C_G_min,connSW,vrb);

            % Price (directed Barabàsi-Albert)
            %fprintf('BA generation\n')
            %n_ini = ceil(rand*n/2);
            n_ini = 5;
            %etaBA = rand;
            etaBA = 0.5;
            m_ini = round(etaBA*n_ini*(n_ini-1));
            %pBA = rand;
            %dnew = 1+floor(pBA*n_ini);
            %dnew = 2;
            %connBA = round(rand);
            connBA = 1;
            BA1 = Price(n_ini,m_ini,n,2,connBA,vrb);
            BA2 = Price(n_ini,m_ini,n,3,connBA,vrb); 
            BA3 = Price(n_ini,m_ini,n,4,connBA,vrb); 

            % saving graphs and their parameters
            % if W == 1
            %     save("graphs_ER_SW_BA.mat","ER2","params_ER2","SW","params_SW","BA","params_BA");
            % end

            % Clustering
            for kk = 1:clusters
                k = pp.ks(kk);
                fprintf(['n = ' num2str(n) '; test = ' num2str(t) ...
                   '; k = ' num2str(k) '\n'])
                Phi_Gs_ER1(t,kk) = netcompl(ER1,k,vrb,type,dic,dic_info);
                Phi_Gs_ER2(t,kk) = netcompl(ER2,k,vrb,type,dic,dic_info);
                Phi_Gs_ER3(t,kk) = netcompl(ER3,k,vrb,type,dic,dic_info);
                Phi_Gs_SW1(t,kk) = netcompl(SW1,k,vrb,type,dic,dic_info);
                Phi_Gs_SW2(t,kk) = netcompl(SW2,k,vrb,type,dic,dic_info);
                Phi_Gs_SW3(t,kk) = netcompl(SW3,k,vrb,type,dic,dic_info);
                Phi_Gs_BA1(t,kk) = netcompl(BA1,k,vrb,type,dic,dic_info);
                Phi_Gs_BA2(t,kk) = netcompl(BA2,k,vrb,type,dic,dic_info);
                Phi_Gs_BA3(t,kk) = netcompl(BA3,k,vrb,type,dic,dic_info);
            end

            
        end
        toc_ = toc;

        fprintf(['This test time: ' num2str(toc_) ' s\n'])
        avg_time_per_test = avg_time_per_test + toc_;
        time_left = round(avg_time_per_test*(tests/t-1)/60);
        if mod(t,1) == 0 && time_left >= 1
            fprintf([num2str(time_left) ' mins left.\n\n'])
        else
            fprintf('\n')
        end
    end
    avg_time_per_test = avg_time_per_test / (60*tests);
    fprintf(['Avg test time: ' num2str(avg_time_per_test)...
        ' mins/test\n\n'])
    

    %% averaging
    % fprintf('Averaging\n')
    % pp.Phi_Gs_ER1 = zeros(1,nodes);
    % pp.Phi_Gs_ER2 = zeros(1,nodes);
    % pp.Phi_Gs_SW = zeros(1,nodes);
    % pp.Phi_Gs_BA = zeros(1,nodes);
    % for n_ = 1:length(pp.N)
    %     pp.Phi_Gs_ER1(n_) = mean(Phi_Gs_ER1(:,n_));
    %     pp.Phi_Gs_ER2(n_) = mean(Phi_Gs_ER2(:,n_));
    %     pp.Phi_Gs_SW(n_) = mean(Phi_Gs_SW(:,n_));
    %     pp.Phi_Gs_BA(n_) = mean(Phi_Gs_BA(:,n_));
    % end

    %% plotting figures
    fprintf('Drawing\n')
    %makesFigures(pp,k,type,tests)
    pp.ftsz = 25;
    makesBoxplots(pp,k,Phi_Gs_ER1,Phi_Gs_ER2,Phi_Gs_ER3,...
        Phi_Gs_SW1,Phi_Gs_SW2,Phi_Gs_SW3,...
        Phi_Gs_BA1,Phi_Gs_BA2,Phi_Gs_BA3);

    fprintf('Done\n')
%end



%% final figure maker
function [] = makesFigures(pp,k,type,tests)

N = pp.N;
lw = pp.lw;
ftsz = pp.ftsz;
inter = pp.inter;
marksz = pp.marksz;

figure('Position',[200 200 600 600])

hER1 = plot(N,pp.Phi_Gs_ER1,'k','linewidth',lw);
hold on
grid on
plot(N,pp.Phi_Gs_ER1,'ok','MarkerSize',marksz,...
    'MarkerEdgeColor','k','MarkerFaceColor','k');

hER2 = plot(N,pp.Phi_Gs_ER2,'r','linewidth',lw);
plot(N,pp.Phi_Gs_ER2,'or','MarkerSize',marksz,...
    'MarkerEdgeColor','r','MarkerFaceColor','r');

hSW = plot(N,pp.Phi_Gs_SW,'g','linewidth',lw);
plot(N,pp.Phi_Gs_SW,'og','MarkerSize',marksz,...
    'MarkerEdgeColor','g','MarkerFaceColor','g');

hBA = plot(N,pp.Phi_Gs_BA,'b','linewidth',lw);
plot(N,pp.Phi_Gs_BA,'ob','MarkerSize',marksz,...
    'MarkerEdgeColor','b','MarkerFaceColor','b');

n = max(N);
xtcks = [1 50:50:200];
xlabel('$n$','Interpreter',inter,'FontSize',ftsz);
xticks(xtcks)
xtickangle(0)
xlim([1 n])
ylabel('$\Phi_{G}(n)$','Interpreter',inter,'FontSize',ftsz);
ymin = 0;
ytcks = 0:50:200;
yticks(ytcks)
ylim([ymin n])
legend([hER1 hER2 hSW hBA], {strcat('Erdos-Renyi $p=0.25 \ln(n)/n',... % ,num2str(pp.pER1)
    '$'),strcat('Erdos-Renyi $p=\ln(n)/n','$'),... % ,num2str(pp.pER2)
    'Small-World', 'Scale-free'},...
    'LineWidth',lw,'FontSize',0.93*ftsz,'Interpreter',inter,...
    'location','northwest');
set(gca,'TickLabelInterpreter',inter,'fontsize',ftsz)
tests_str = '$ tests.';
if tests == 1
    tests_str = '$ test.';
end
title(['CL: ' type ', $k = ' num2str(k) '$, averaged over $'...
    num2str(tests) tests_str],...
    'interpreter',inter,'FontSize',ftsz/1.25);
drawnow

end




%% final figure maker
function [] = makesBoxplots(pp,k,Phi_Gs_ER1,Phi_Gs_ER2,Phi_Gs_ER3,...
    Phi_Gs_SW1,Phi_Gs_SW2,Phi_Gs_SW3,...
    Phi_Gs_BA1,Phi_Gs_BA2,Phi_Gs_BA3)

N = pp.N;
lw = pp.lw;
ftsz = pp.ftsz;
inter = pp.inter;
marksz = pp.marksz;

%figure('Position',[200 200 600 600])
%boxplot([Phi_Gs_ER1, Phi_Gs_ER2, Phi_Gs_ER3, Phi_Gs_SW1, Phi_Gs_SW2, Phi_Gs_SW3, Phi_Gs_BA1, Phi_Gs_BA2, Phi_Gs_BA3],...
    %{'ER thr = 0.25 ln(n)/n','ER thr = 0.5 ln(n)/n','ER thr = 1.25 ln(n)/n',...
    %'SW p = 0.25','SW p = 0.5','SW p = 0.75',...
    %'SF dnew = 2','SF dnew = 3','SF dnew = 4'})


figure('Position',[100 100 1100 400])
boxplot([Phi_Gs_SW1(:,1), Phi_Gs_SW1(:,2), Phi_Gs_SW2(:,1), Phi_Gs_SW2(:,2), Phi_Gs_SW3(:,1), Phi_Gs_SW3(:,2)],...
    {'$p = 0.25$, $k = 2$','$p = 0.25$, $k = 4$','$p = 0.5$, $k = 2$','$p = 0.5$, $k = 4$','$p = 0.75$, $k = 2$','$p = 0.75$, $k = 4$'})
set(gca,'TickLabelInterpreter',inter,'fontsize',ftsz)
yticks(85:5:100)
ylabel('$\varphi_{\mathcal{G}}$','Interpreter',inter,'FontSize',ftsz);

%title([' n = ' num2str(N)]) 

figure('Position',[200 100 1100 400])
boxplot([Phi_Gs_BA1(:,1), Phi_Gs_BA1(:,2), Phi_Gs_BA2(:,1), Phi_Gs_BA2(:,2), Phi_Gs_BA3(:,1), Phi_Gs_BA3(:,2)],...
    {'$m = 2$, $k = 2$','$m = 2$, $k = 4$','$m = 3$, $k = 2$','$m = 3$, $k = 4$','$m = 4$, $k = 2$','$m = 4$, $k = 4$'})
set(gca,'TickLabelInterpreter',inter,'fontsize',ftsz)
yticks(50:10:100)
ylabel('$\varphi_{\mathcal{G}}$','Interpreter',inter,'FontSize',ftsz);


%title([' n = ' num2str(N)]) % 'k = ' num2str(k) 
drawnow

end



