    function [a1, a2] = twoAxes(parent)
        gl = uigridlayout(parent, [2 1]);
        a1 = uiaxes(gl);
        a2 = uiaxes(gl);
    end
