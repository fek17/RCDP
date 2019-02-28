%% extents of reaction
figure
yyaxis left
hold on
for i = 1:size(c.RX,1)
    t.this = sprintf('xi_%u',i);
    plot(v.z,v.(t.this),'DisplayName',t.this)
end
ylabel('\xi_i / kmol.h^{-1}');
xlabel('z / m');
yyaxis right
plot(v.z,t.y(:,6),'DisplayName','T / K')
ylabel('T / K');
legend('Location','east');

% export
formatFig(12,12);
% print(gcf, '-dpdf', [pwd '/graphs/overview-extents.pdf']);

%% molar flows

figure
yyaxis left
hold on
for i = 1:numel(c.species)
    plot(v.z,v.(['n_' c.species{i}]),'DisplayName',['n_{' c.species{i} '}'])
end
ylabel('n_j / kmol.h^{-1}');
xlabel('z / m');
yyaxis right
plot(v.z,t.y(:,6),'DisplayName','T / K')
ylabel('T / K');
legend('Location','east');

% export
formatFig(12,12);
% print(gcf, '-dpdf', [pwd '/graphs/overview-molar-flows.pdf']);

%% yield of PA wrt OX

figure
plot(v.z,v.Y_PA_OX)
ylabel('Y_{PA,OX}^{OV}');
xlabel('z / m');
