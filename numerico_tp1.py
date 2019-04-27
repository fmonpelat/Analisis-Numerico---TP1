#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
TP1 ANALISIS NUMERICO
Rocío Gallo
Facundo Monpelat
Problemas de busqueda de raices
Comparacion de metodos: Newton Raphson, Biseccion y Punto Fijo
m0 = 100000/97490 = 1.025
Lo = 2 * 100000/97490 = 2.051
err = 0.5*10ˆ-15 
"""

from scipy.optimize import brentq
import numpy as np #Manejo de arrays
import sys
sys.setrecursionlimit(10000000)
import math 
# para graficar
import plotly
import plotly.plotly as py
plotly.tools.set_credentials_file(username='fmonpelat', api_key='YpD7z4O34340q0GZGbC7')
#plotly.tools.set_credentials_file(username='fmonpelat', api_key='xxxxxx')
import plotly.graph_objs as go
# debug
from pprint import pprint

# funcion lineal
def fl(y): 
    return y
def dfl(y): 
    return 1
def ddfl(y): 
    return 0

# Funcion TP cuando m=0 (ayuda con las derivadas con el wolfram alpha)
def f(y):
    return (-2*10*y*(1-(2.051)/math.sqrt(math.pow(y,2)+1)))
def gf(y):
    return (y*2.051/math.sqrt(math.pow(y,2)+1)) 
    #return (y-2*10*y*(1-(2.051)/math.sqrt(math.pow(y,2)+1)))
def df(y): 
    return (-20-((41.02*math.pow(y,2))/math.pow((1 + math.pow(y,2)),3/2)) + 41.02/math.sqrt(1+math.pow(y,2)))
def ddf(y): 
    return ((-123.06*y)/((math.sqrt(1 + math.pow(y,2) ) *(1. + math.pow(y,2) )**2)))

# Funcion TP cuando tenemos mg (ayuda con las derivadas con el wolfram alpha)
#def f2(y,m=1.025*1.5): 
def f2(y,m=0): 
    g=9.81
    k=10
    Lo=2.051
    a=1
    #return (-2*k*y-m*g+((2*k*y*Lo)/(math.sqrt(math.pow(y,2)+math.pow(a,2)))))
    return (-2*10*y - m*g + (( 2*10*y*2.051)/math.sqrt(math.pow(y,2)+1)))
    #(-2*10*y - 0.3075*9.81 + (( 2*10*y*2.051)/sqrt(pow(y,2)+1)))

# primera funcion de g(x) para el metodo de punto fijo ej2
def gf2(y,m=1.025*0.3):
    g=9.81
    k=10
    Lo=2.051
    a=1
    #return (y + 2*k*y - (2*k*y*Lo)/(math.sqrt(math.pow(y,2)+math.pow(a,2))) + m*g)
    #return ( (m*g)/(-2*k+(2*k*Lo)/(math.sqrt(math.pow(y,2)+math.pow(a,2)))) )

    # este es el return por la cual las 2 raices de los extremos funciona.
    #return ((y*Lo/(math.sqrt(math.pow(y,2)+math.pow(a,2)))) - (m*g)/(2*k))
    #(1.025*9.81)/((-2*10*2.051)/sqrt(y**2+1))
    #return ((20*y +9.81*1.025*0.3)*(math.sqrt(y**2+1)))/(2*10*2.051) 
    # este return nos da la ec de g para la raiz positiva primera (dar el int donde esta g funciona)
    return ((-m*g)/(2*k-2*k*Lo/(math.sqrt(math.pow(y,2)+math.pow(a,1) ))))

# segunda funcion de g(x) para el metodo de punto fijo ej2
def gf3(y,m=1.025*0.3):
    g=9.81
    k=10
    Lo=2.051
    a=1
    #return (y + 2*k*y - (2*k*y*Lo)/(math.sqrt(math.pow(y,2)+math.pow(a,2))) + m*g)
    #return ( (m*g)/(-2*k+(2*k*Lo)/(math.sqrt(math.pow(y,2)+math.pow(a,2)))) )
    # este es el return por la cual las 2 raices de los extremos funciona.
    return ((y*Lo/(math.sqrt(math.pow(y,2)+math.pow(a,2)))) - (m*g)/(2*k))
    

def df2(y): 
    #con los parametros:
    #g=9.81
    #k=10
    #Lo=2.051
    #a=1
    #m0=1.025 (IMPORTANTE: la derivada es la misma independientenente de su masa (su termino cuando se deriva es 0))
    return (-20-((41.02*math.pow(y,2))/(math.pow((1 + math.pow(y,2)),3/2))) + 41.02/math.sqrt(1 + math.pow(y,2)))
    #return (-20-(41.02*math.pow(y,2))/(1 + math.pow(math.pow(y,2),3/2)) + 41.02/math.sqrt(1 + math.pow(y,2)) )




#Funcion que imprime gráfico mediante plotly de las diferencias entre muestras
# Recibe:
# NewtonRaphsonData = datos en array de cada diferencia por el metodo de newton
# fixedPointData = datos en array de cada diferencia por el metodo de punto fijo
# BisectData = datos en array de cada diferencia por el metodo de bisección
def GraficarDiferencia(NewtonRaphsonData,fixedPointData,BisectData):
    # data1,data2,data3 son arrays con puntos a ser graficados
    # TODO hay que ajustar estos valores de linspace ya que el gráfico es logaritmico
    xx = np.linspace(0, 100, 100)

    trace0 = go.Scatter(
                        x=xx,
                        y=fixedPointData,
                        name = 'Datos del método Punto Fijo',
                        mode = 'lines+markers',
    )
    trace1 = go.Scatter(
                        x = xx,
                        y = NewtonRaphsonData,
                        name = 'Datos del método Newton Raphson',
                        mode = 'lines+markers',
    )
    trace2 = go.Scatter(
                    x = xx,
                    y = BisectData,
                    name = 'Datos del método de Bisección',
                    mode = 'lines+markers',
    )
    layout = dict(title = 'Gráfica de diferencias entre muestras por metodo',
                xaxis = dict(title = 'X (iteraciones)'),
                yaxis = dict(title = 'Y (diferencia entre muestras)',
                             type='log',
                             autorange=True),
                )
    data = [trace0,trace1,trace2]
    #la figura es un dicccionario de data(traces) y el layout (que tambien es un dicccionario)
    fig = dict(data=data, layout=layout)
    plotly.offline.plot(fig, auto_open=True)
    return

#Funcion que devuelve un intervalo de convergencia para alguna raiz de una funcion con el metodo de newton raphson
# Recibe:
# f = Función
# df = Serivada de la función
# nRoot = Semilla inicial (se utiliza como semilla la raiz desde donde nos movemos)
# intPrecision = Incremento en el eje x 
# counter_array = Array devuelto por referencia que contendra las iteraciones maximas hasta obtener el intervalo.
def intConv(f,df,nRoot,intPrecision,counter_array):
    root = nRoot
    pos = 0
    n    = 0
    start = 0
    end = 0
    end_i = 0
    start_i = 0
    maxXpos = 50
    n_max=1000 # si hay 1000 iteraciones y no se llego al valor de precision a_tol se toma como que el metodo diverge
    a_tol = 0.5e-15
    arr_delta = []
    printdata = []
    debug = 0
    
    pos = nRoot
    while( abs(root-nRoot) < 0.1 ):
      pos = pos + intPrecision
      root,n = NewtonRaphson(f,df,pos,a_tol,n,n_max,arr_delta,printdata)
      if (debug==1): print("cicle: {0:.2f} iteration: {1:.2f} posx: {2:.2f} root: {3:.2f}".format(end_i,n,pos,root))
      end_i=end_i+1
      n=0
      end = pos
      if (debug==1): print("dif root:{0:.2f} pos:{1:.2f} nRoot:{2:.2f} dif:{3:.2f}".format(root,pos,nRoot,abs(root-nRoot)))
      if(abs(nRoot-pos)>maxXpos):
        break
    if (debug==1): print("posx: {0:.2f}".format(pos))

    pos = nRoot
    root = nRoot
    while( abs(root-nRoot) < 0.1 ):
      pos = pos - intPrecision
      root,n = NewtonRaphson(f,df,pos,a_tol,n,n_max,arr_delta,printdata)
      if (debug==1): print("cicle: {0:.2f} iteration: {1:.2f} posx: {2:.2f} root: {3:.2f}".format(start_i,n,pos,root))
      start_i=start_i+1
      n=0
      start = pos
      if(abs(nRoot-pos)>maxXpos):
        break
    counter_array.append(start_i)
    counter_array.append(end_i)
    return start, end
    
    """
    # iteracion lado positivo
    while (math.fabs(nEnd-nStart) > intPrecision):
        n=0
        nMiddle=(nEnd-nStart)/2

        print("middle:{0:5f}".format(nMiddle))
        try:
            root,n = NewtonRaphson(f,df,nMiddle,a_tol,n,n_max,arr_delta,printdata)
        except ValueError:
            print("divirgió! middle: {0:.5f}".format(nMiddle))
            # divirgio la funcion devolvio la exepcion value error (movemos a la izquierda el puntero de end)
            nEnd = nMiddle
        if (math.fabs(root-nRoot) < 0.01): 
            # converge a la misma raiz pasada con un error de comparacion de 0.01 (movemos el puntero de start a la derecha)
            nStart=nMiddle
        else:
            print("convergio a otra raiz! middle: {0:.5f} root: {1:.5f} nRoot: {2:.5f}".format(nMiddle,root,nRoot))
            # convergió a otra raiz (debemos mover el puntero de start en la raiz que converge y la raiz suministrada dividido 2)
            nEnd=nRoot
            nStart = root
            print("new nEnd: {0:.5f} new nStart: {1:.5f}".format(nEnd,nStart))
    nStInt = nStart
    print("Start :{0:.10f}".format(nStInt))
    return nStInt, nEnInt
    """

#Funcion que imprime el intervalo de convergencia mediante el metodo de newton
# Recibe:
# f2 = funcion 
# df2 = derivada de la funcion f2
# root = raiz para el intervalo
def PrintMaxInterval(f2,df2,root):
    IncPrecision = 0.0001
    maxXpos=50
    start, end = intConv(f2,df2,root,IncPrecision,counter_array)
    if(abs(start)>maxXpos):
      print("Intervalo de maxima convergencia: [start,end]: Error+/-{2:.4f} [{0:},{1:.4f}]".format("infinity",end,IncPrecision))
    elif(abs(end)>maxXpos):
      print("Intervalo de maxima convergencia: [start,end]: Error+/-{2:.4f} [{0:.4f},{1:}]".format(start,"infinity",IncPrecision))
    else:
      print("Intervalo de maxima convergencia: [start,end]: Error+/-{2:.4f} [{0:.4f},{1:.4f}]".format(start,end,IncPrecision))

#Funciones busqueda de raices
def bisec(function, start, stop, a_tol, n_max, arr_delta,printdata):
    """
    Devolver (x0, delta), raiz y cota de error por metodo de la biseccion
    Datos deben cumplir f(a)*f(b) > 0
    """
    debug = 0
    i = 0
    x = start+(stop-start)/2    #mejor que (a+b)/2 segun Burden
    delta = (stop-start)/2
    #print('{0:^4} {1:^17} {2:^17} {3:^17}'.format('i', 'x', 'start', 'stop'))
    #print('{0:4} {1: .14f} {2: .14f} {3: .14f}'.format(0, x, start, stop))
    
    for i in range(n_max):
        if function(start) * function(x) > 0:
            start = x
        else:
            stop = x
        x_old = x
        x = start+(stop-start)/2 #(a+b)/2
        delta = np.abs(x - x_old)
        #print('{0:4} {1: .14f} {2: .14f} {3: .14f}'.format(i+1, x, start, stop))
        aux_dict={
        'k':i,
        'ak':start,
        'bk':stop,
        'Fak':function(start),
        'Fbk':function(stop),
        'Rk1':x,
        'deltaRk1':delta,
        'deltaR/R':delta/x,
        }
        printdata.append(aux_dict)
        arr_delta.append(delta)
        if delta <= a_tol: #Hubo convergencia
            if (debug==1): print('Hubo convergencia, n_iter = ' + str(i+1))
            return x, delta, i+1

    #Si todavia no salio es que no hubo convergencia:
    raise ValueError('No hubo convergencia')
    return x, delta, i+1

def FixedPoint(f,seed,a_tol,n,n_max,arr_delta,printdata):
    debug = 1
    #  para que converga la seed debe satisfacer |f'(seed)|<1
    # arrdelta devuelve un array de diferencias entre valores
    if n>n_max: 
        raise ValueError('No hubo convergencia')
        return root_aux, n

    root_aux = f(seed)
    delta = np.abs(root_aux - seed)
    # la impresion se delega devolviendo los valores como array de diccionarios y usando la función de impresión.
    #print('K,Yk,delta Y,deltaY/Yk,lambda,p')
    #print ('{0:.1f},{1:.15f},{2:.15f},{3:.15f},{4:},{5:}'.format(n,root_aux,delta,delta/root_aux,"N/a","N/a"))
    aux_dict={
        'k':n,
        'root':root_aux,
        'delta':delta,
        'delta/root':delta/root_aux,
    }
    printdata.append(aux_dict)
    arr_delta.append(delta)
    if delta <= a_tol:
        if (debug==1): print('Hubo convergencia, n_iter = ' + str(n))
        return root_aux, n
    else:
        return FixedPoint(f,root_aux,a_tol,n+1,n_max,arr_delta,printdata)

def NewtonRaphson(f,df,seed,a_tol,n,n_max,arr_delta,printdata):
    debug = 1
    #  para que converga la seed debe satisfacer |f'(seed)|<1
    if n>n_max: 
        raise ValueError('No hubo convergencia')
        return root_aux, n

    root_aux = seed - f(seed)/df(seed)

    delta = np.abs(root_aux - seed)
    #print('K,Yk,delta Y,deltaY/Yk,lambda,p')
    #print ('{0:.1f},{1:.4f},{2:.4f},{3:.4f},{4:.4f}'.format(n,root_aux,delta,delta/root_aux,seed))
    aux_dict={
        'k':n,
        'root':root_aux,
        'delta':delta,
        'delta/root':delta/root_aux,
    }
    printdata.append(aux_dict)
    if delta <= a_tol:
        arr_delta.append(delta)
        if (debug==1): print('Hubo convergencia, n_iter = ' + str(n))
        return root_aux, n
    else:
        arr_delta.append(delta)
        return NewtonRaphson(f,df,root_aux,a_tol,n+1,n_max,arr_delta,printdata)

def secante(f, x0, x1, a_tol, n_max):
    """
    Devolver (x, delta), raiz y cota de error por metodo de la secante
    """
    delta = 0
    debug = 0
    print('{0:^4} {1:^17} {2:^17} {3:^17}'.format('i', 'x', 'x_-1', 'delta'))
    print('{0:4} {1: .14f} {2: .14f} {3: .14f}'.format(0, x1, x0, delta))

    for i in range(n_max):
        x = x1 - f(x1)*(x1-x0)/(f(x1)-f(x0))
        delta = np.abs(x - x1)
        x0 = x1
        x1 = x
        
        print('{0:4} {1: .14f} {2: .14f} {3: .14f}'.format(i+1, x1, x0, delta))
        
        #Chequear convergencia
        if delta <= a_tol: #Hubo convergencia
            if (debug==1): print('Hubo convergencia, n_iter = ' + str(i+1))
            return x, delta, i+1

    #Si todavia no salio es que no hubo convergencia:
    raise ValueError('No hubo convergencia')
    return x. delta, i+1

#funciones de impresión
def print_convergence_table(printdata):
    ## truncamiento a primeros 10 iteraciones y ultimas 10 iteraciones
    trunc_iter=10
    print ('---------------------------------------------------')
    print('k, $Y_{k}$, $\Delta Y$, $\dfrac{\Delta Y}{Y_{k}}$')
    for i in range(len(printdata)):
      k=printdata[i]['k']+1
      print ('{0:},{1:.15f},{2:.15f},{3:.15f}'.format(k,printdata[i]['root'],printdata[i]['delta'],printdata[i]['delta/root']))

    ## latex style
    print ('--- LATEX TABLE  ---')
    #printing only hline rows and 10 upmost and 10 leftmost
    length = len(printdata)
    if( length+1  > trunc_iter ):
      start=trunc_iter
      stop= length - trunc_iter
      del printdata[start:stop]
    #printing only hline rows
    print('\\begin{table}[H]')
    print('\makegapedcells\n\centering\n\\resizebox{0.55\\textwidth}{!}{')
    print('\\begin{tabular}{|c|c|c|c|}')
    print('\hline k & Y$_{k}$ & $\Delta$ Y & $\dfrac{\Delta Y}{Y_{k}}$ \\\\')
    for i in range(len(printdata)):
      k=printdata[i]['k']+1
      print ('\hline {0:} & {1:.15f} & {2:.15f} & {3:.15f} \\\\'
      .format(k,printdata[i]['root'],printdata[i]['delta'],printdata[i]['delta/root']))
    print('\hline\n\end{tabular}}')
    print('\end{table}')
    print ('---------------------------------------------------')
    return

def print_start_table(printdata):
    ## truncamiento a primeros 10 iteraciones y ultimas 10 iteraciones
    trunc_iter=10
    print ('---------------------------------------------------')
    print('k, ak, bk, f(ak), f(bk), rk+1, deltar k+1, delta r/r')
    for i in range(len(printdata)):
      k=printdata[i]['k']+1
      print ('{0:},{1:.4f},{2:},{3:.4f},{4:.4f},{5:.15f},{6:.4f},{7:.4f}'
      .format(k,printdata[i]['ak'],printdata[i]['bk'],printdata[i]['Fak'],printdata[i]['Fbk'],printdata[i]['Rk1'],printdata[i]['deltaRk1'],
      printdata[i]['deltaR/R']))

    ## latex style
    print ('--- LATEX TABLE  ---')
    #printing only hline rows and 10 upmost and 10 leftmost
    length = len(printdata)
    if( length+1  > trunc_iter ):
      start=trunc_iter
      stop= length - trunc_iter
      del printdata[start:stop]

    print('\\begin{table}[H]')
    print('\makegapedcells\n\centering\n\\resizebox{0.55\\textwidth}{!}{')
    print('\\begin{tabular}{|c|c|c|c|c|c|c|c|}')
    print('\hline k & a$_{k}$ & b$_{k}$ & f(a$_{k}$) & f(b$_{k}$) & r$_{k+1}$ & $\Delta r_{k+1}$ & $\dfrac{\Delta r}{r}$ \\\\')
    for i in range(len(printdata)):
      k=printdata[i]['k']+1
      print ('\hline {0:} & {1:.4f} & {2:.4f} & {3:.4f} & {4:.4f} & {5:.16f} & {6:.4f} & {7:.4f} \\\\'
      .format(k,printdata[i]['ak'],printdata[i]['bk'],printdata[i]['Fak'],printdata[i]['Fbk'],printdata[i]['Rk1'],printdata[i]['deltaRk1'],
      printdata[i]['deltaR/R']))
    print('\hline\n\end{tabular}}')
    print('\end{table}')
    print ('---------------------------------------------------')
    return



##########
## EJ1
#Encontrar los puntos de equilibrio del sistema sin tener en cuenta el efecto de la
#gravedad. Para esto alcanza con considerar m=0 en la ecuación no lineal. Hallar el punto de
#equilibrio positivo mediante los siguientes métodos con error absoluto menor a 0.5 10-15
#(en módulo).
##########
print ('--------------------------------------------------------------------------------')
print('--- EJ 1 ---')

a_tol1 = 0.5e-15
n_max = 100

print('----------------')
print('Metodo biseccion')
print('----------------')
#Intervalo para buscar raiz
a = 0.2
b = 3
data3 = [] # array para guardar los deltas del metodo en cada iteracion
printdata=[]
r, delta, n_iter = bisec(f, a, b, a_tol1, n_max,data3,printdata)
print_start_table(printdata)


print('----------------')
print('Metodo punto fijo')
print('----------------')
iterator=0
root=0.72
data2 = []
printdata=[]
print('K,Yk,delta Y,deltaY/Yk')
root, iterator = FixedPoint(gf,root,a_tol1,iterator,n_max,data2,printdata)
print_convergence_table(printdata)


print('----------------')
print('Metodo Newton Raphson')
print('----------------')
iterator=0
root=-1
data1 = []
printdata=[]
print('K  ,Yk,delta Y,deltaY/Yk')
root, iterator = NewtonRaphson(f,df,root,a_tol1,iterator,n_max,data1,printdata)
print_convergence_table(printdata)

GraficarDiferencia(data1,data2,data3)



##########
## EJ2
#Repetir el procedimiento realizado en el ítem 1) para el caso m=0.3*m0
#Encontrar todos los puntos de equilibrio del sistema. Presentar tablas y gráficos análogos al ítem 1) y
#comparar resultados. Tomar como criterio de corte del método la diferencia absoluta entre
#dos iteraciones consecutivas menor a 0.5x10-15
##########
print ('--------------------------------------------------------------------------------')
print('--- EJ 2 ---')

a_tol1 = 0.5e-15
n_max = 10000

print('----------------')
print('Metodo biseccion')
print('----------------')
#Intervalo para buscar raiz
a = -3
b = -1
dataB1 = [] # array para guardar los deltas del metodo en cada iteracion
dataB2 = []
dataB3 = []
printdata=[]

r, delta, n_iter = bisec(f2, a, b, a_tol1, n_max,dataB1,printdata)
print_start_table(printdata)
print ('-----------------')
a = -0.5
b = 0.6
printdata=[]
r, delta, n_iter = bisec(f2, a, b, a_tol1, n_max,dataB2,printdata)
print_start_table(printdata)
print ('-----------------')
a = 0.8
b = 5
printdata=[]
r, delta, n_iter = bisec(f2, a, b, a_tol1, n_max,dataB3,printdata)
print_start_table(printdata)

print('----------------')
print('Metodo punto fijo')
print('----------------')
iterator=0
dataPF1 = []
dataPF2 = []
dataPF3 = []
printdata=[]

#print('K,Yk,delta Y,deltaY/Yk,lambda,p')
root=-1
root, iterator = FixedPoint(gf3,root,a_tol1,iterator,n_max,dataPF1,printdata)
print_convergence_table(printdata)
print ('-----------------')

root=0
iterator=0
printdata=[]
# cambiamos de funcion g(x) ya que no converge a la raiz primera negativa en el intervalo proximo a la raiz
root, iterator = FixedPoint(gf2,root,a_tol1,iterator,n_max,dataPF2,printdata)
print_convergence_table(printdata)
print ('-----------------')
root=1
iterator=0
printdata=[]
root, iterator = FixedPoint(gf3,root,a_tol1,iterator,n_max,dataPF3,printdata)
print_convergence_table(printdata)


print('----------------')
print('Metodo Newton Raphson')
print('----------------')
iterator = 0
dataNR1 = []
dataNR2 = []
dataNR3 = []
printdata = []
counter_array = []
#print('K, Yk, delta Y, deltaY/Yk, lambda, p')

root = -2
printdata=[]
root, iterator = NewtonRaphson(f2,df2,root,a_tol1,iterator,n_max,dataNR1,printdata)
#PrintMaxInterval(f2,df2,root)
print_convergence_table(printdata)
print ('-----------------')
root = 0.33
iterator = 0
root, iterator = NewtonRaphson(f2,df2,root,a_tol1,iterator,n_max,dataNR2,printdata)
#PrintMaxInterval(f2,df2,root)
print_convergence_table(printdata)
print ('-----------------')
root = 1
iterator = 0
printdata=[]
root, iterator = NewtonRaphson(f2,df2,root,a_tol1,iterator,n_max,dataNR3,printdata)
#PrintMaxInterval(f2,df2,root)
print_convergence_table(printdata)

GraficarDiferencia(dataNR1,dataPF1,dataB1)
GraficarDiferencia(dataNR2,dataPF2,dataB2)
GraficarDiferencia(dataNR3,dataPF3,dataB3)



##########
## EJ3
#Solamente para el método de Newton-Raphson y aplicado al caso del ítem 2), encontrar el
#máximo intervalo de convergencia de cada raíz, es decir, el intervalo de mayor longitud
#posible en el cual cualquier semilla perteneciente a él hace que el método converja
##########
print ('--------------------------------------------------------------------------------')
print('--- EJ 3 ---')

iterator = 0
dataNR1 = []
dataNR3 = []
counter_array = []
root = -1
printdata=[]
root, iterator = NewtonRaphson(f2,df2,root,a_tol1,iterator,n_max,dataNR1,printdata)
PrintMaxInterval(f2,df2,root)

root = 1
iterator = 0
printdata=[]
root, iterator = NewtonRaphson(f2,df2,root,a_tol1,iterator,n_max,dataNR3,printdata)
PrintMaxInterval(f2,df2,root)



