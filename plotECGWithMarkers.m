function ax = plotECGWithMarkers(ecg, Fs, titleStr, varargin)
%PLOTECGWITHMARKERS Plot ECG with multiple marker sets + optional time range.
%
%   plotECGWithMarkers(ecg, Fs, titleStr, ...
%       "tmin", 0, "tmax", 3.5, ...
%       idx1, "label1", "rv", idx2, "label2", "ko", ...)
%
% Marker sets can be:
%   idx, label
%   idx, label, style   (e.g. "rv", "k*", "go")
%
% Name–Value options (can be anywhere before markers):
%   "tmin", "tmax", "Parent"

arguments
    ecg {mustBeVector, mustBeNumeric}
    Fs (1,1) double {mustBePositive, mustBeFinite}
    titleStr {mustBeTextScalar}
end

arguments (Repeating)
    varargin
end

ecg = ecg(:);
N   = numel(ecg);
time = (0:N-1)'/Fs;
titleStr = string(titleStr);

% ---- Parse leading name–value options (tmin/tmax/Parent) ----
tmin = [];
tmax = [];
ax   = [];

knownNames = ["tmin","tmax","parent"];
k = 1;

while k <= numel(varargin)
    if (ischar(varargin{k}) || isstring(varargin{k})) && any(lower(string(varargin{k})) == knownNames)
        name = lower(string(varargin{k}));
        if k+1 > numel(varargin)
            error("Option '%s' has no value.", name);
        end

        switch name
            case "tmin"
                tmin = varargin{k+1};
            case "tmax"
                tmax = varargin{k+1};
            case "parent"
                ax = varargin{k+1};
        end

        k = k + 2;
    else
        break
    end
end

markerArgs = varargin(k:end);

% ---- Create axes if needed ----
if isempty(ax) || ~isgraphics(ax,"axes")
    figure;
    ax = axes;
end

hold(ax,"on");
grid(ax,"on");

% ---- Plot ECG ----
plot(ax, time, ecg, "DisplayName","ECG");

% ---- Parse marker sets: (idx,label[,style]) ----
specs = struct("idx",{}, "label",{}, "style",{});
i = 1;
while i <= numel(markerArgs)
    if i+1 > numel(markerArgs)
        error("Markers must be passed as idx,label[,style] sets.");
    end

    idx   = markerArgs{i};
    label = markerArgs{i+1};

    style = ""; % default: auto
    if i+2 <= numel(markerArgs) && (ischar(markerArgs{i+2}) || isstring(markerArgs{i+2})) ...
            && isLineSpec(markerArgs{i+2})
        style = string(markerArgs{i+2});
        i = i + 3;
    else
        i = i + 2;
    end

    specs(end+1) = struct("idx",idx, "label",string(label), "style",style); %#ok<AGROW>
end

% ---- Plot markers ----
colors       = lines(10);
markerStyles = ["o","s","^","v","d",">","<","p","h","x","+"];

mCount = 0;

for s = 1:numel(specs)
    idx = specs(s).idx;
    if isempty(idx), continue; end

    idx = idx(:);
    idx = double(idx);
    valid = isfinite(idx) & idx >= 1 & idx <= N;
    idx = unique(round(idx(valid)));
    if isempty(idx), continue; end

    mCount = mCount + 1;
    col = colors(mod(mCount-1, size(colors,1)) + 1, :);
    mk  = markerStyles(mod(mCount-1, numel(markerStyles)) + 1);

    if specs(s).style ~= ""
        % User provided style -> use it
        h = plot(ax, time(idx), ecg(idx), specs(s).style, ...
            "LineStyle","none", ...
            "MarkerSize",6, ...
            "DisplayName", specs(s).label);

        % If marker is fillable, fill it with its own color
        if isFillableMarker(specs(s).style)
            set(h, "MarkerFaceColor","auto");
        end
    else
        % Auto style
        plot(ax, time(idx), ecg(idx), mk, ...
            "LineStyle","none", ...
            "MarkerSize",6, ...
            "MarkerFaceColor", col, ...
            "Color", col, ...
            "DisplayName", specs(s).label);
    end
end

xlabel(ax,"Time [s]");
ylabel(ax,"ECG Amplitude [mV]");
title(ax, titleStr);
legend(ax,"show","Location","best");

% ---- Apply x-limits ----
xl = xlim(ax);
if ~isempty(tmin), xl(1) = tmin; end
if ~isempty(tmax), xl(2) = tmax; end
xlim(ax, xl);

end


% ===================== local helpers =====================

function tf = isLineSpec(s)
% Heuristic: typical short MATLAB linespecs like "rv", "k*", "go", "b."
s = string(s);
tf = (strlength(s) <= 6) && ~contains(s," ");
if ~tf, return; end

% Must contain at least one marker or line char or color char
markers = ".ox+*sdv^<>ph";
lines   = "-:";
colors  = "rgbcmykw";
tf = any(contains(s, split(markers,""))) || any(contains(s, split(lines,""))) || any(contains(s, split(colors,"")));
end

function tf = isFillableMarker(style)
% Markers that can be filled: o s d ^ v > < p h
style = char(style);
tf = any(ismember(style, ['o','s','d','^','v','>','<','p','h']));
end
