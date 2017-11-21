figure(1); hold on;
k = 3; 
for V = 0:0.05:5
    polycoeff = [-1 (V-k) 0 0 -1 V];
    rts = roots(polycoeff);
    rts = rts(imag(rts) == 0);
    plot(V*ones(length(rts),1),rts,'r.');
end
hold off;
xlabel('V'); ylabel('Zeros points');