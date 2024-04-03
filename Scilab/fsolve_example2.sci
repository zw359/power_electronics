// Example 2 
func1 = 'res(1)=x(2)-(x(1).^2+1)';
func2 = 'res(2)=x(1)-(2*x(2)-x(2).^2)/3';
//deff('res=fct_2(x)',['res(1)=x(2)-(x(1).^2+1)';'res(2)=x(1)-(2*x(2)-x(2).^2)/3']) 
deff('res=fct_2(x)',[func1; func2]) 
scf(2) 
clf(2) 
x1 = linspace(-3,3,101) 
y1 = x1.^2+1 
y2 = linspace(-3,5,51) 
x2=(2*y2-y2.^2)/3 
plot(x1,y1,'r-'); 
p = get("hdl"); p.children.thickness = 3; 
plot(x2,y2,'b-'); 
p = get("hdl"); p.children.thickness = 3; 

x0 = [1;0]
xsol1 =fsolve(x0,fct_2) 
res1 = fct_2(xsol1) 

x0 = [-3;8]
xsol2 =fsolve(x0,fct_2) 
res2 = fct_2(xsol2)
