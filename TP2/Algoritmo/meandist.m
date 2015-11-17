%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%Funcion brindada por la catedra de Metodos Numericos Avanzados, ITBA
%Parametros de entrada:
% * data = dataset
% * code = code vectors de un individuo
%
%Parametros de salida:
% * dist = distorsi ́on media
function dist = meandist(data,code)
[dim,n] = size(data);
[dim,nvector] = size(code);
dist = 0;
for l = 1:n
%busco la mınima distancia/distorsion de cada columna
    distances = norm(repmat(data(:,l),1,nvector)-code,2,'columns');
    [m,i] = min(distances);
    dist = dist+m;
end
dist = dist/n;