classdef TestClass < handle
    properties
        Var1;
        Var2;
    end
    
    methods
        function obj = TestClass(var1init, var2init)
            obj.Var1 = var1init;
            obj.Var2 = var2init;
        end
        
        function [var1, var2] = update(obj)
            obj.Var1 = obj.Var1*2;
            obj.Var2 = obj.Var2*2;
            var1 = obj.Var1;
            var2 = obj.Var2;
        end
    end
end