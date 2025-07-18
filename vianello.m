clc;
clear;
clear all;

n = input('Matrisin boyutunu girin (n x n için n değeri): ');

fprintf('\nMatris elemanlarını satır satır giriniz:\n');
A = zeros(n);  
for i = 1:n
    for j = 1:n
        A(i,j) = input(sprintf('A(%d,%d) = ', i, j));
    end
end

fprintf('\nBaşlangıç vektörünü giriniz:\n');
x0 = zeros(n,1);  
for i = 1:n
    x0(i) = input(sprintf('x0(%d) = ', i));
end

tolerans = input('\nHata toleransını girin : ');
maxiterasyon = 1000;  

%% En Büyük Özdeğeri Bulma (Vianello Metodu)
lambda_buyuk_eski = 0;
hata = inf;
iterasyon = 0;
x = x0/norm(x0);  

while (hata > tolerans) && (iterasyon < maxiterasyon)
    iterasyon = iterasyon + 1;
    y = A * x; 
    lambda_buyuk_yeni = norm(y, inf);  
    x = y / lambda_buyuk_yeni;  
    hata = abs(lambda_buyuk_yeni - lambda_buyuk_eski);
    lambda_buyuk_eski = lambda_buyuk_yeni;
end

fprintf('\nEn büyük özdeğer (lambda_max): %.8f\n', lambda_buyuk_yeni);

%% En Küçük Özdeğeri Bulma (Matris Tersiyle)
if det(A) == 0
    fprintf('\nUYARI: Matrisin determinantı sıfır, en küçük özdeğer bulunamaz!\n');
else
    A_tersi = inv(A); 
    lambda_kucuk_eski = 0;
    hata = inf;
    iterasyon = 0;
    x = x0/norm(x0);  
    
    while (hata > tolerans) && (iterasyon < maxiterasyon)
        iterasyon = iterasyon + 1;
        y = A_tersi * x;
        lambda_tersi_yeni = norm(y, inf);
        x = y / lambda_tersi_yeni;
        hata = abs(lambda_tersi_yeni - lambda_kucuk_eski);
        lambda_kucuk_eski = lambda_tersi_yeni;
    end
    lambda_kucuk_yeni = 1 / lambda_tersi_yeni;  
    fprintf('En küçük özdeğer (lambda_min): %.8f\n', lambda_kucuk_yeni);
end