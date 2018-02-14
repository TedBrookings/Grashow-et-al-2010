function TestMapLocations
numMaps = 30;
k = 3;

mapList = generateMapList(numMaps, k);
clustInds = kmeans_pp(mapList, k, 10, @mapDistance)

return

function mapList = generateMapList(numMaps, k)
protoMapList = [LocationMap(round(rand(7,7)))];
for n = 2:k
  protoMapList = [protoMapList ; LocationMap(round(rand(7,7)))];
end


mapList = [protoMapList(1)];
for n = 2:numMaps
  mapList = [mapList ; protoMapList(1 + mod(n-1,k))];
end
return