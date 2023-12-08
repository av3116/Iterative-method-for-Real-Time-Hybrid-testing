function bH7 = copyBroyden(bH1)
bH7 = eval(class(bH1));  %create default object of the same class as a. one valid use of eval
p_all =  properties(bH1).';
for i=1:length(p_all)  %copy all public properties
    p = p_all{i};
    try   %may fail if property is read-only
        bH7.(p) = bH1.(p);
    catch
        warning('failed to copy property: %s', p);
    end
end
end