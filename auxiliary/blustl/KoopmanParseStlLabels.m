function stl = KoopmanParseStlLabels(Sys)
% STLC_parse_stl_labels     Parses the STL specifications of an STLC_lti
%                           instance and replaces state, output, input and
%                           disturbance labels with the corresponding indices
%                           in X, Y, U and W.
%
% Input:
%       Sys: an STLC_lti instance
%
% Output:
%       stlList: STL specification over X, Y, U and W.
%
% :copyright: TBD
% :license: TBD

labels = {'xlabel','ulabel'};
vars = {'x','u'};

stl=Sys.stl;
for ilabel = 1:numel(labels)
    var = vars{ilabel};
    for jlabel = 1:numel(Sys.(labels{ilabel}))
        label = Sys.(labels{ilabel}){jlabel};
        stl = regexprep(stl,['\<' label '\(t\)'], [ var '(' num2str(jlabel) ',t)'] );
    end
end
