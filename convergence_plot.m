function convergence_plot(history, x0, x, par)
% SLP iteration history:

iters = 1:length(history);
cons_viol = arrayfun(@(h) h.constraint_violation, history);
duality_gap = arrayfun(@(h) h.duality_gap, history);
delta_vals = arrayfun(@(h) h.delta, history);
step = arrayfun(@(h) h.step, history);
F_comp_vals = arrayfun(@(h) h.F_comp_val, history);
F_surp_vals = arrayfun(@(h) h.F_surp_val, history);
quality = arrayfun(@(h) h.quality, history);

font = 26;

figure()
subplot(3,1,1)
semilogy(iters, cons_viol, '-+k','LineWidth',2.5, 'MarkerSize', 12);
ylabel('Residual','Interpreter','latex');
xlim([1,length(iters)])
title('SLP Convergence','Interpreter','latex');
set(gca,'Fontsize',font)

subplot(3,1,2)
semilogy(iters, abs(duality_gap), '-+k', 'LineWidth', 2.5, 'MarkerSize', 12);
ylabel('Duality','Interpreter','latex');
xlim([1,length(iters)])
set(gca,'Fontsize',font)

subplot(3,1,3)
plot(iters, delta_vals, '-+k', 'LineWidth', 2.5, 'MarkerSize', 12);
hold on
semilogy(iters, step, '--+r', 'LineWidth', 2.5, 'MarkerSize', 12);
xlabel('Iteration','Interpreter','latex');
ylabel('Step','Interpreter','latex');
xlim([1,length(iters)])
lgd = legend({'$\|\overline{\delta x}\|_{\infty}$', '$\|\delta x\|_{\infty}$'}, 'Interpreter', 'latex', 'Orientation', 'horizontal');
set(lgd, 'Color', 'none');
set(lgd, 'Box', 'off');
set(lgd, 'Units', 'pixels');
set(lgd, 'Position', [-27,10.323,400,33.353]);
set(gca, 'FontSize', 24);
set(gca,'Fontsize',font)
set(gcf,'units','points','position',[18,240,607,625])

figure()
subplot(3,1,1)
plot(iters, F_comp_vals * par.varq, '-+k', 'LineWidth', 2.5, 'MarkerSize', 12);
ylabel('Energy','Interpreter','latex');
title('Objective Values','Interpreter','latex');
xlim([1,length(iters)])
set(gca,'Fontsize',font)

subplot(3,1,2)
plot(iters, F_surp_vals, '-+k', 'LineWidth', 2.5, 'MarkerSize', 12);
xlabel('Iteration','Interpreter','latex');
ylabel('Supply','Interpreter','latex');
xlim([1,length(iters)])
set(gca,'Fontsize',font)

subplot(3,1,3)
plot(iters, quality, '-+k', 'LineWidth', 2.5, 'MarkerSize', 12);
xlabel('Iteration','Interpreter','latex');
ylabel('Demand','Interpreter','latex');
xlim([1,length(iters)])
set(gca,'Fontsize',font)
set(gcf,'units','points','position',[626,240,605,625])


nP = length(par.index.nodes);
nvar = length(x0);

idxP = 1:nP;
idxF = nP+1:nvar;

end