function varargout = LoadAbf(FileName, varargin)
%AbfStruct = LoadAbf(FileName, WantedString)

[t,data,h] = load_abf(FileName);
NumChan = h.nADCNumChannels;
if(size(data, 2) ~= NumChan)
  data = permute(data, [2 1 3]);
end

OutStruct.Time = t;
OutStruct.Header = h;
NumFound = 0;
for m = 1:NumChan
  SampleNum = h.nADCSamplingSeq(m) + 1;
  FieldName = deblank(h.sADCChannelName(SampleNum,:));
  Units = deblank(h.sADCUnits(SampleNum,:));
  FieldData = squeeze(data(:,m,:));

  if(nargin > 1)
    NotWanted = true;
    for n = 1:length(varargin)
      if(length(strfind(FieldName, varargin{n})) > 0)
	NotWanted = false;
	break
      end
    end
    if(NotWanted)
      continue
    end
  end
  
  NumFound = NumFound + 1;
  OutStruct.Data.(FieldName) = FieldData;
  OutStruct.Units.(FieldName) = Units;
end

if(NumFound == 0)
  error('Unable to find requested wave data')
end

if(nargout == 1)
  varargout = {OutStruct};
  return
end

%When is the rest of this shit ever used????

varargout = {t};

FNames = fieldnames(OutStruct.Data);
for n = 1:length(FNames)
  if(strcmp(OutStruct.Units.(FNames{n}), 'mV'))
    if(nargin > 1)
      for m=1:length(varargin)
	if(length(strfind(FNames{n}, varargin{m})) > 0)
	  varargout = {varargout{:}, OutStruct.Data.(FNames{n})};
	  break;
	end
      end
    else
      varargout = {varargout{:}, OutStruct.Data.(FNames{n})};
    end
  end
end

if(nargout > length(varargout))
  for n = 1:length(FNames)
    if(strcmp(OutStruct.Units.(FNames{n}), 'nA'))
      if(nargin > 1)
	for m=1:length(varargin)
	  if(length(strfind(FNames{n}, varargin{m})) > 0)
	    varargout = {varargout{:}, OutStruct.Data.(FNames{n})};
	    break;
	  end
	end
      else
	varargout = {varargout{:}, OutStruct.Data.(FNames{n})};
      end
    end
  end
end