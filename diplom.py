from sympy import *
x = symbols('x')
init_printing(use_unicode=True)


def residual(a,b,func,result):
    t=a
    h=(b-a)/100
    max_residual=0
    while(t <= b):        
        temp=abs(func.subs(x,t)-result.subs(x,t))
        if(temp>max_residual):
            max_residual=temp
        t+=h
    return max_residual
def main (start,end,mu,eps_mu,func):
    z_a = start
    z_b = end
    temp_a=start
    temp_b=end
    while(True):
        
        f0=func.subs(x,z_a)
        f1=func.subs(x,z_b)

        d0=dfunc.subs(x,z_a)
        d1=dfunc.subs(x,z_b)

        alpha=(d0*(z_b**3 - z_a**3)-3*(f1-f0)*z_a**2)/((z_b-z_a)*(5*z_a**3+5*z_b*z_a**2 +2*z_a*z_b**2))
        beta=(d1*(z_b**3 - z_a**3)-3*(f1-f0)*z_b**2)/((z_b-z_a)*(5*z_b**3+5*z_a*z_b**2 +2*z_b*z_a**2))

        lamda=(2*z_a**2-z_b**2-z_b*z_a)/(5*z_a**3+5*z_b*z_a**2 +2*z_a*z_b**2)
        gamma=(2*z_b**2-z_a**2-z_b*z_a)/(5*z_b**3+5*z_a*z_b**2 +2*z_b*z_a**2)

        c=(alpha-beta)/(gamma-lamda)
        b=alpha+c*lamda
        a=(f1-f0+b*(z_a**2-z_b**2)+c*(z_a-z_b))/(z_b**3-z_a**3)
        d=f0-a*z_a**3-b*z_a**2-c*z_a
            
        result = a*x**3 +b*x**2+c*x+d
        max_mu=residual(z_a,z_b,func,result)
        print("max = {0}, a = {1}, b = {2}".format(max_mu,z_a,z_b))
        
        if (max_mu > mu):
            z_b = (temp_a + temp_b)/ 2
            temp_b=z_b
           
            continue    
        elif ((max_mu < mu*(1-eps_mu))&(temp_b < x1)&(max_mu!=0)):
            temp = z_b
            z_b = (3*temp_b-temp_a)/2
            temp_b=z_b
            temp_a=temp       
            continue
        else:
            if(z_b == x1):
    #             print("max = {0}, a = {1}, b = {2}".format(max_mu,z_a,z_b))
                print("")
                break        
    #         print("max = {0}, a = {1}, b = {2}".format(max_mu,z_a,z_b))
            print("----------------------------------------")
            z_a = z_b
            z_b = x1
            temp_a = z_a
            temp_b = z_b
            continue
    
    




# Write Your function
# func = exp(x)
func=sin(x)
x0=0.5
x1=3.1

#Set residual!!
mu=0.02
eps_mu=0.03
#Out
print(" Вигляд функції f(x):") 
pprint(func)

print("mu={0}, eps_mu={1}, a={2}, b={3}".format(mu,eps_mu,x0,x1))

print("\nЗагальний вигляд сплайна F(A,x):")
A,B,C,D = symbols('A B C D')
pprint(A*x**3 +B*x**2+C*x+D)
print("\n")
 
dfunc=diff(func)


# f0=func.subs(x,x0)
# f1=func.subs(x,x1)

# d0=dfunc.subs(x,x0)
# d1=dfunc.subs(x,x1)

# alpha=(d0*(x1**3 - x0**3)-3*(f1-f0)*x0**2)/((x1-x0)*(5*x0**3+5*x1*x0**2 +2*x0*x1**2))
# beta=(d1*(x1**3 - x0**3)-3*(f1-f0)*x1**2)/((x1-x0)*(5*x1**3+5*x0*x1**2 +2*x1*x0**2))

# lamda=(2*x0**2-x1**2-x1*x0)/(5*x0**3+5*x1*x0**2 +2*x0*x1**2)
# gamma=(2*x1**2-x0**2-x1*x0)/(5*x1**3+5*x0*x1**2 +2*x1*x0**2)

# c=(alpha-beta)/(gamma-lamda)
# b=alpha+c*lamda
# a=(f1-f0+b*(x0**2-x1**2)+c*(x0-x1))/(x1**3-x0**3)
# d=f0-a*x0**3-b*x0**2-c*x0



z_a = x0
z_b = x1
temp_a=x0
temp_b=x1
while(True):
    
    f0=func.subs(x,z_a)
    f1=func.subs(x,z_b)

    d0=dfunc.subs(x,z_a)
    d1=dfunc.subs(x,z_b)

    alpha=(d0*(z_b**3 - z_a**3)-3*(f1-f0)*z_a**2)/((z_b-z_a)*(5*z_a**3+5*z_b*z_a**2 +2*z_a*z_b**2))
    beta=(d1*(z_b**3 - z_a**3)-3*(f1-f0)*z_b**2)/((z_b-z_a)*(5*z_b**3+5*z_a*z_b**2 +2*z_b*z_a**2))

    lamda=(2*z_a**2-z_b**2-z_b*z_a)/(5*z_a**3+5*z_b*z_a**2 +2*z_a*z_b**2)
    gamma=(2*z_b**2-z_a**2-z_b*z_a)/(5*z_b**3+5*z_a*z_b**2 +2*z_b*z_a**2)

    c=(alpha-beta)/(gamma-lamda)
    b=alpha+c*lamda
    a=(f1-f0+b*(z_a**2-z_b**2)+c*(z_a-z_b))/(z_b**3-z_a**3)
    d=f0-a*z_a**3-b*z_a**2-c*z_a
        
    result = a*x**3 +b*x**2+c*x+d
    max_mu=residual(z_a,z_b,func,result)
    print("max = {0}, a = {1}, b = {2}".format(max_mu,z_a,z_b))
    
    if (max_mu > mu):
        z_b = (temp_a + temp_b)/ 2
        temp_b=z_b
       
        continue    
    elif ((max_mu < mu*(1-eps_mu))&(temp_b < x1)&(max_mu!=0)):
        temp = z_b
        z_b = (3*temp_b-temp_a)/2
        temp_b=z_b
        temp_a=temp       
        continue
    else:
        if(z_b == x1):
#             print("max = {0}, a = {1}, b = {2}".format(max_mu,z_a,z_b))
            print("")
            break        
#         print("max = {0}, a = {1}, b = {2}".format(max_mu,z_a,z_b))
        print("----------------------------------------")
        z_a = z_b
        z_b = x1
        temp_a = z_a
        temp_b = z_b
        continue
    
# # print("max ={:4f}").format(max_mu)
# print("a=",z_a,"b=",z_b)
# #--------------------------------
# # Out

# print(" Вигляд функції f(x):") 
# pprint(func)

# print("\nЗагальний вигляд сплайна F(A,x):")
# A,B,C,D = symbols('A B C D')
# pprint(A*x**3 +B*x**2+C*x+D)

# print("\n Результати обчислень")
# print("A ={0:0.5f}".format(a))
# print("B ={0:0.5f}".format(b))
# print("C ={0:0.5f}".format(c))
# print("D ={0:0.5f}".format(d))

# result = a*x**3 +b*x**2+c*x+d

# print("\n Похибка апроксимації")
# max_residual = residual(x0,x1,func,result)
# print("max ={0:0.5f}".format(max_residual))
# #---------------------------------
# p1=plot(func,(x,x0*0.98,z_b*1.02),show=False,line_color='b',axis_center="auto",
#         title=" Графік функції та сплайна",ylabel="")
# p2=plot(result,(x,x0*0.98,x1*1.02),show=False,line_color='r')
# p1.extend(p2)
# p1.show()