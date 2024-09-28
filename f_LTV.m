function F=f_LTV(aprime,hprime,a,h,z)

F = 0;

if aprime<0 && hprime~=0
    F = -aprime/hprime;
end

end %end function "Chen2010_ReturnFn"
