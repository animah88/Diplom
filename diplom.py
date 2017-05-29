from sympy import *
import numpy as np
x = symbols('x')
init_printing(use_unicode = True)

def toString(arr):
    return list(map(lambda x: str(x), arr))

def myRound(after_point):
    def f(val):
        return round(val, after_point)
    return f

def residual(a,b,func,result):
    t = a
    h = (b-a)/51
    max_residual = 0

    while(t <= b):        
        temp = abs(func.subs(x,t)-result.subs(x,t))
        # print("SHIT!",result.subs(x,t) )
        if(temp > max_residual):

            max_residual=temp
        t+=h
    return max_residual

pres = myRound(3)

def main (start,end,mu,eps_mu,func,typ):
    print("Type is:",typ)
    z_a = start
    z_b = end

    x_val_func = np.linspace(start,end,50)
    y_val_func = np.vectorize(lambdify(x,func))(x_val_func)

    temp_a=start
    temp_b=end

    dfunc=diff(func)
    array = []
    while(True):
        
        f0=func.subs(x,z_a)
        f1=func.subs(x,z_b)

        d0=dfunc.subs(x,z_a)
        d1=dfunc.subs(x,z_b)

        if typ == "polynomial":            
            alpha = (d0*(z_b**3 - z_a**3)-3*(f1-f0)*z_a**2)/((z_b-z_a)*(5*z_a**3+5*z_b*z_a**2 +2*z_a*z_b**2))
            beta = (d1*(z_b**3 - z_a**3)-3*(f1-f0)*z_b**2)/((z_b-z_a)*(5*z_b**3+5*z_a*z_b**2 +2*z_b*z_a**2))

            lamda = (2*z_a**2-z_b**2-z_b*z_a)/(5*z_a**3+5*z_b*z_a**2 +2*z_a*z_b**2)
            gamma = (2*z_b**2-z_a**2-z_b*z_a)/(5*z_b**3+5*z_a*z_b**2 +2*z_b*z_a**2)

            c=N((alpha-beta)/(gamma-lamda))
            b=N(alpha+c*lamda)
            a=N((f1-f0+b*(z_a**2-z_b**2)+c*(z_a-z_b))/(z_b**3-z_a**3))
            d=N(f0-a*z_a**3-b*z_a**2-c*z_a)
                
            result = a*x**3 +b*x**2+c*x+d

        if typ == "exponential":
            A = ( (d0*z_a*ln(z_b/z_a)-f0*ln(f1/f0)) / (f0*(2*ln(z_b/z_a)*z_a**2 - z_b**2 +z_a**2)) )/((f1*(z_a*ln(z_b/z_a)-z_b+z_a)*(2*ln(z_b/z_a)*z_b**2-z_b**2+z_a**2)) - (( z_a*ln(z_b/z_a)-z_b+z_a)*(2*ln(z_b/z_a)*z_b**2-z_b**2+z_a**2))/(z_b*ln(z_b/z_a)-z_b+z_a )*(2*ln(z_b/z_a)*z_a**2-z_b**2+z_a**2))
            B = (z_b*d1*ln(z_b/z_a)-ln(f1/f0)*f1) / (f1*(2*ln(z_b/z_a)*z_b**2-z_b**2+z_a**2)*(f1*(2*ln(z_b/z_a)*z_a**2-z_b**2+z_a**2)*(z_b*ln(z_b/z_a)-z_b+z_a)-1))
            
            c = N(A-B)
            b = N(((z_b*d1*ln(z_b/z_a))/(f1*(z_b*ln(z_b/z_a)-z_b+z_a))) - c*((2*ln(z_b/z_a)*z_b**2-z_b**2+z_a**2)/(z_b*ln(z_b/z_a)-z_b+z_a))-ln(f1/f0)/(z_b*ln(z_b/z_a)-z_b+z_a))
            d = N((ln(f1/f0)-b*(z_b-z_a)-c*(z_b**2-z_a**2))/(ln(z_b/z_a)))
            a = N(f0/(exp(b*z_a+c*z_a**2)*z_a**d))
            print("SHIT!",a,b,c,d )
            result = a*exp(b*x+c*x**2)*(x**d)


        max_mu = residual(z_a,z_b,func,result)
        # print("max = {0}, a = {1}, b = {2}".format(max_mu,z_a,z_b))
        
        if (max_mu > mu):
            z_b = (temp_a + temp_b)/ 2
            temp_b=z_b
           
            continue    
        elif ((max_mu < mu*(1-eps_mu))&(temp_b < end)&(max_mu!=0)):
            print("max = {0}, a = {1}, b = {2}".format(max_mu,z_a,z_b))
            temp = z_b
            z_b = (3*temp_b-temp_a)/2
            temp_b=z_b
            temp_a=temp       
            continue
        else:
            # array.append([z_a,z_b,a,b,c,d])
            t=z_a
            h=(z_b-z_a)/50
            array_x = []
            array_y = []
            residual_x=[]
            residual_y=[]
            while(t <= z_b):     
                array_x.append(t)
                array_y.append(result.subs(x,t))
                residual_x.append(t)
                # residual_y.append(abs(result.subs(x,t)-func.subs(x,t)))
                residual_y.append(abs(result.subs(x,t)-func.subs(x,t)))
                # print("residual",abs(result.subs(x,t)-func.subs(x,t)))
                # temp=abs(func.subs(x,t)-result.subs(x,t))
                t+=h
            array_x.append(z_b)
            residual_x.append(z_b)
            array_y.append(result.subs(x,z_b))
            residual_y.append(abs(result.subs(x,z_b)-func.subs(x,z_b)))
             
            if (typ=="polynomial"):
                latech = latex(pres(a)*x**3 + pres(b)*x**2 + pres(c)*x + pres(d))
            else:
                latech = latex(pres(a)*(exp( pres(b)*x + pres(c)*x**2  ))*(x**pres(d)))

            array.append({
                'start': pres(z_a),
                'end': pres(z_b),
                'formula': latech,
                'x': list(map(lambda x: float(x), array_x)),
                'y': list(map(lambda y: float(y), array_y)),
                'func_x': list(x_val_func),
                'func_y': list(y_val_func),
                'residual_x': list(map(lambda x: float(x), residual_x)),
                'residual_y': list(map(lambda y: float(y), residual_y))
                })

            if(z_b == end):
                print("max = {0}, a = {1}, b = {2}".format(max_mu,z_a,z_b))
                print("")
                break        
            print("max = {0}, a = {1}, b = {2}".format(max_mu,z_a,z_b))
            print("----------------------------------------")
            z_a = z_b
            z_b = end
            temp_a = z_a
            temp_b = z_b
            continue


    
    # return list(map(toString,array))
    return array
    
    

