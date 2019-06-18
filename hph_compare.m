% Compares the model covariance (the first column of HPH') with the data
% covariance (xcov(X)) to see if we are, in fact, using a realistic model.
% Creates a new plot every time.

function hph_compare(H,P,X)
HPH = hph(H,P);
xcovX = xcov(X,'unbiased');
xcovX = flip(fftshift(xcovX));
xcovX = xcovX(1:((length(xcovX)+1)/2));
figure
plot(HPH,'.-');hold on
plot(xcovX,'.-')
legend('HPH^T','covariance')
xlabel('lag')
end
