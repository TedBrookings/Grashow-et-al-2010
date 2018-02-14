classdef LocationMap
   properties (SetAccess = public)
     mapMatrix
   end
   %properties (SetAccess = protected)
   %  %
   %end
   methods
     %constructor
     function obj = LocationMap(mapMatrix)
       obj.mapMatrix = mapMatrix;
     end
     %operators
     function mapC = plus(mapA, mapB)
       mapC = LocationMap(mapA.mapMatrix + mapB.mapMatrix);
     end
     function mapC = minus(mapA, mapB)
       mapC = LocationMap(abs(mapA.mapMatrix - mapB.mapMatrix));
     end
     function divMap = rdivide(map, scalar)
       divMap = LocationMap(round(map.mapMatrix / scalar));
     end
     function divMap = mrdivide(map, scalar)
       divMap = LocationMap(round(map.mapMatrix / scalar));
     end     

     %overloaded functions
     function fbool = isnane(mapList)
       fbool = repmat(false, size(mapList));
     end
     function fbool = isfinite(mapList)
       fbool = repmat(true, size(mapList));
     end
     function sumMap = sum(mapList)
       sumMap = mapList(1);
       for n = 2:length(mapList)
	 sumMap = sumMap + mapList(n);
       end
     end
     function meanMap = mean(mapList)
       meanMap = sum(mapList) / length(mapList);
       meanMap.mapMatrix = round(meanMap.mapMatrix);
     end
     function mapMat = repmat(map, numRows, numCols)
       mapCol = map;
       for n = 2:numCols
	 mapCol = [mapCol , map];
       end
       mapMat = mapCol;
       for n = 2:numRows
	 mapMat = [mapMat ; mapCol];
       end
     end
     
     %class-specific functions
     function d = mapDistance(mapA, mapB)
       d = 0;
       [burstRowsA, burstColsA] = find(mapA.mapMatrix == 1);
       [burstRowsB, burstColsB] = find(mapB.mapMatrix == 1);
       for nA = 1:length(burstRowsA)
	 %d = d + min(abs(burstRowsB - burstRowsA(nA)) + ...
		%     abs(burstColsB - burstColsA(nA)));
	 d = d + min(sqrt((burstRowsB - burstRowsA(nA)).^2 + ...
			  (burstColsB - burstColsA(nA)).^2));
       end
       for nB = 1:length(burstRowsB)
	 %d = d + min(abs(burstRowsA - burstRowsB(nB)) + ...
		%     abs(burstColsA - burstColsB(nB)));
	 d = d + min(sqrt((burstRowsA - burstRowsB(nB)).^2 + ...
			  (burstColsA - burstColsB(nB)).^2));
       end       
       d = d / (length(burstRowsA) + length(burstRowsB));
     end
   end
end



