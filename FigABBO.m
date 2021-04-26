%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Codes for estimating lattice constant of ABB'O using PLSR and PCR.
% Dataset can be seen in data_ABBO.mat. Dataset contains experimental
% lattice constant [C1], ionic radii [C2-C5], atomic number [C6-C9],
% electronegativity [C10-C13], and density [C14-C17].
%
% Output:
% estimated constant: yfitReg (OLSR), yfitPLS (PLSR), yfitPCR (PCR)
% estimation error: error_reg (OLSR), error_plsr (PLSR), error_pcr (PCR)
% variance explained: PctVar (PLSR), PctExpl (PCR)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load Data
load('data_ABBO.mat')
y_train = data(:, 1);
X_train = data(:, 2:end);

%% Plot parameters
% comp2plot = number of components to consider for regression
n1 = size(X_train, 1);
fsize = 20;
pcpair = 4;
comp2plot = [4,8,11];
num2plot = size(comp2plot, 2);
f1 = figure( 'Units', 'normalized', 'Position', [0.1 0.25 0.75 1] );

%% OLS Regression
b = regress(y_train, X_train);
yfitReg = X_train*b;
error_reg = mean(abs((y_train - yfitReg) ./ y_train)*100);

%% Print and store results
fprintf('Regression = %0.4f\n', error_reg)
fprintf('No. of Comp \t PLS \t \t \t PCR\n')

for i = comp2plot
    %% PLSR
    [XL,~,~,~,betaPLS, PctVar, MSE] = plsregress(X_train,y_train,i);
    yfitPLS = [ones(n1,1) X_train]*betaPLS;
    
    %% PCR
    [PCALoadings,PCAScores,~,~,PctExpl] = pca(X_train);
    betaPCR = regress(y_train-mean(y_train), PCAScores(:,1:i));
    betaPCR = PCALoadings(:,1:i)*betaPCR;
    betaPCR = [mean(y_train) - mean(X_train)*betaPCR; betaPCR];
    yfitPCR = [ones(n1,1) X_train]*betaPCR;
    
    %% Estimation error and cumulative variance explained
    error_plsr = mean(abs((y_train - yfitPLS) ./ y_train)*100);
    error_pcr = mean(abs((y_train - yfitPCR) ./ y_train)*100);
    
    cumvar_plsr = 100*sum(PctVar(1,1:i));
    cumvar_pcr = sum(PctExpl(1:i));
    
    %% Print and plot results
    fprintf('%d \t \t PLS = %0.3f \t PCR = %0.3f\n',...
        i, error_plsr, error_pcr)
    
    figure(f1)
    subplot(2,num2plot, round(i/pcpair))
    plot(y_train, yfitPLS, 'o', 'MarkerSize', 15);
    hold on
    eb_l = plot(7.6:10, 7.6:10, '--', 'Color', [0.6 0.6 0.6],...
        'LineWidth', 5);
    title({['PLSR: ', num2str(i), ' Comp'],...
        sprintf('Error: %0.4f%%' , error_plsr),...
        sprintf('CumVar: %0.4f%%' , cumvar_plsr)})
    xlabel('Experimental Value')
    ylabel('Estimated Value')
    ylim([7.6 9.2])
    set(gca, 'ytick', 7.6:0.4:9.2, 'FontSize', fsize)
    
    subplot(2,num2plot, round(i/pcpair)+num2plot)
    plot(y_train, yfitPCR, 'o', 'MarkerSize', 15);
    hold on
    eb_l = plot(7.6:10, 7.6:10, '--', 'Color', [0.6 0.6 0.6],...
        'LineWidth', 5);
    title({['PCR: ', num2str(i), ' Comp'],...
        sprintf('Error: %0.4f%%' , error_pcr),...
        sprintf('CumVar: %0.4f%%' , cumvar_pcr)})
    xlabel('Experimental Value')
    ylabel('Estimated Value')
    ylim([7.6 9.2])
    set(gca, 'ytick', 7.6:0.4:9.2, 'FontSize', fsize)
    
    %% Plot loadings for regression with i=8 components
    if i == 8
        f2 = figure( 'Units', 'normalized', 'Position', [0.1 0.25 1 1] );
        for pci = 1:i
            subplot(1,i,pci)
            imagesc(abs(XL(:, pci)))
            set(gca, 'ytick', 1:16,...
                'yticklabel', {'r_A', 'r_B', 'r_{B^\prime}', 'r_O',...
                '\chi_A', '\chi_B', '\chi_{B^\prime}', '\chi_O',...
                'Z_A', 'Z_B', 'Z_{B^\prime}', 'Z_O',...
                '\rho_A', '\rho_B', '\rho_{B^\prime}', '\rho_O'},...
                'xtick', 1, 'xticklabel', ['Comp ', num2str(pci)],...
                'FontSize', fsize)
            title({'[PLSR]', 'Var Exp:',...
                sprintf('%0.2f%%' , PctVar(1, pci)*100)})
            
            % print loading values in plot
            textStrings = num2str(XL(:, pci), '%0.2f');
            textStrings = strtrim(cellstr(textStrings));
            [xplot, yplot] = meshgrid(1, 1:16);
            hStrings = text(xplot(:), yplot(:), textStrings(:),...
                'HorizontalAlignment', 'center', 'FontSize', fsize);
        end
        
        f3 = figure( 'Units', 'normalized', 'Position', [0.1 0.25 1 1] );
        for pci = 1:i
            subplot(1,i,pci)
            imagesc(abs(PCALoadings(:, pci)))
            set(gca, 'ytick', 1:16,...
                'yticklabel', {'r_A', 'r_B', 'r_{B^\prime}', 'r_O',...
                '\chi_A', '\chi_B', '\chi_{B^\prime}', '\chi_O',...
                'Z_A', 'Z_B', 'Z_{B^\prime}', 'Z_O',...
                '\rho_A', '\rho_B', '\rho_{B^\prime}', '\rho_O'},...
                'xtick', 1, 'xticklabel', ['Comp ', num2str(pci)],...
                'FontSize', fsize)
            title({'[PCR]', 'Var Exp:', sprintf('%0.2f%%' , PctExpl(pci))})
            
            % print loading values in plot
            textStrings = num2str(PCALoadings(:, pci), '%0.2f');
            textStrings = strtrim(cellstr(textStrings));
            [xplot, yplot] = meshgrid(1, 1:16);
            hStrings = text(xplot(:), yplot(:), textStrings(:),...
                'HorizontalAlignment', 'center', 'FontSize', fsize);
        end
    end
end