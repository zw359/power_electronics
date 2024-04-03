// Example 1
deff('res=fct_1(x)','res=cos(x)-x^2-x')
x0 = -2;

xsol =fsolve(x0,fct_1)
x = linspace(-2,2,51)
fcos = cos(x)
fx = x.*x+x
scf(1)
clf(1)
plot(x,fcos,'r-');
p = get("hdl"); p.children.thickness = 3;
plot(x,fx,'b-');
p = get("hdl"); p.children.thickness = 3;
