global NN NN
global U UI 
figure(2)
hold on

t = linspace(0,400*(7/400),400);
for pp =1:length(U)
    plot(t,(U(pp,:)/0.001)/1,'k-')
    hold on
    text(1,1,num2str(pp))
    pause(0.05)
    hold off
end