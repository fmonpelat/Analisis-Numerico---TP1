#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
TP1 ANALISIS NUMERICO
Rocío Gallo
Facundo Monpelat
Problemas de busqueda de raices
Comparacion de metodos
m0 = 100000/97490 = 1.025
Lo = 2 * 100000/97490 = 2.051
err = 0.5*10ˆ-15 
"""

#Imports
from scipy.optimize import brentq
import timeit #Para calcular tiempo de corrida
import numpy as np #Manejo de arrays
import sys
sys.setrecursionlimit(500)
# para graficar
import plotly
import plotly.plotly as py
plotly.tools.set_credentials_file(username='fmonpelat', api_key='YpD7z4O34340q0GZGbC7')
import plotly.graph_objs as go


#Funciones prueba
def f1(x): 
    return x**2-2
def df1(x): 
    return 2*x
def ddf1(x): 
    return 2
def gf1(x):
    return (x+2)/(x+1)

# funcion lineal
def fl(x): 
    return x
def dfl(x): 
    return 1
def ddfl(x): 
    return 0

# Funcion TP 
def f(t): 
    return 
def df(t): 
    return 0
def ddf(t): 
    return 0

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


def Graficar(root,iterator,trace0descr,trace1descr,graphtitle):
    # para graficar
    xx = np.linspace(0, 4, 256)
    yy1 = f1(xx)
    yy2 = fl(xx)

    trace0 = go.Scatter(
                        x=xx,
                        y=yy1,
                        name = trace0descr,
                        line = dict(
                                    color = ('rgb(205, 12, 24)'),
                                    width = 4
                                    )
                        )
    trace1 = go.Scatter(
                        x = xx,
                        y = yy2,
                        name = trace1descr,
                        line = dict(
                                    color = ('rgb(22, 96, 167)'),
                                    width = 4,
                                    dash = 'dash'
                                    )
                        )
    layout = dict(title = graphtitle,
                xaxis = dict(title = 'X'),
                yaxis = dict(title = 'Y'),
                shapes = [{
                'type': 'line',
                'x0': root,
                'y0': -4,
                'x1': root,
                'y1': 10,
                'line': {
                    'color': 'rgb(55, 128, 191)',
                    'dash': 'dot',
                    'width': 3,}
                }],
                annotations=[
                                dict(
                                    x=root,
                                    y=0,
                                    xref='x',
                                    yref='y',
                                    text='<b>Raiz = '+str(root)+'</b><br>Iteraciones = '+str(iterator),
                                    showarrow=True,
                                    font=dict(
                                        family='Courier New, monospace',
                                        size=15,
                                        color='#000000'
                                    ),
                                    align='center',
                                    arrowhead=2,
                                    arrowsize=1,
                                    arrowwidth=2,
                                    arrowcolor='#636363',
                                    ax=150,
                                    ay=-150,
                                    bordercolor='#ffffff',
                                    borderwidth=2,
                                    borderpad=4,
                                    bgcolor='#ffffff',
                                    opacity=0.8
                                )
        ]
                )
    data = [trace0,trace1]
    #la figura es un dicccionario de data(traces) y el layout (que tambien es un dicccionario)
    fig2 = dict(data=data, layout=layout)
    plotly.offline.plot(fig2, auto_open=True)
    return


def maxIntervalNewtonRaphson(initial,f,df,a_tol):
    # Funcion para el punto 3
    # metodo para encontrar cotas de semillas max donde no diverja el método
    # se toma que diverge cuando n>1000 iteraciones
    #para el lado positivo de la semilla
    n_max=1000
    arr_delta=[]
    n=0
    seed = initial
    root,n = NewtonRaphson(f,df,seed,a_tol,n,n_max,arr_delta)
    return
    

#Funciones busqueda de raices
def bisec(function, start, stop, a_tol, n_max, arr_delta):
    """
    Devolver (x0, delta), raiz y cota de error por metodo de la biseccion
    Datos deben cumplir f(a)*f(b) > 0
    """
    x = start+(stop-start)/2    #mejor que (a+b)/2 segun Burden
    delta = (stop-start)/2
    
    print('{0:^4} {1:^17} {2:^17} {3:^17}'.format('i', 'x', 'start', 'stop'))
    print('{0:4} {1: .14f} {2: .14f} {3: .14f}'.format(0, x, start, stop))
    
    for i in range(n_max):
        if function(start) * function(x) > 0:
            start = x
        else:
            stop = x
        x_old = x
        x = start+(stop-start)/2 #(a+b)/2
        delta = np.abs(x - x_old)
        arr_delta.append(delta)
        print('{0:4} {1: .14f} {2: .14f} {3: .14f}'.format(i+1, x, start, stop))
        
        if delta <= a_tol: #Hubo convergencia
            print('Hubo convergencia, n_iter = ' + str(i+1))
            return x, delta, i+1
    
    #Si todavia no salio es que no hubo convergencia:
    raise ValueError('No hubo convergencia')
    return x, delta, i+1


def FixedPoint(f,seed,a_tol,n,n_max,arr_delta):

    #  para que converga la seed debe satisfacer |f'(seed)|<1
    # arrdelta devuelve un array de diferencias entre valores
    if n>n_max: 
        raise ValueError('No hubo convergencia')
        return root_aux, n

    root_aux = f(seed)
    delta = np.abs(root_aux - seed)
    #print('K,Yk,delta Y,deltaY/Yk,lambda,p')
    print ('{0:.1f},{1:.4f},{2:.4f},{3:.4f},{4:.4f}'.format(n,root_aux,delta,delta/root_aux,seed))
    arr_delta.append(delta)
    if delta <= a_tol:
        print('Hubo convergencia, n_iter = ' + str(n))
        return root_aux, n
    else:
        return FixedPoint(f,root_aux,a_tol,n+1,n_max,arr_delta)


def NewtonRaphson(f,df,seed,a_tol,n,n_max,arr_delta):
    #  para que converga la seed debe satisfacer |f'(seed)|<1
    if n>n_max: 
        raise ValueError('No hubo convergencia')
        return root_aux, n

    root_aux = seed - (f(seed)/df(seed))
    delta = np.abs(root_aux - seed)
    #print('K,Yk,delta Y,deltaY/Yk,lambda,p')
    print ('{0:.1f},{1:.4f},{2:.4f},{3:.4f},{4:.4f}'.format(n,root_aux,delta,delta/root_aux,seed))
    arr_delta.append(delta)
    if delta <= a_tol:
        print('Hubo convergencia, n_iter = ' + str(n))
        return root_aux, n
    else:
        return NewtonRaphson(f,df,root_aux,a_tol,n+1,n_max,arr_delta)



def secante(f, x0, x1, a_tol, n_max):
    """
    Devolver (x, delta), raiz y cota de error por metodo de la secante
    """
    delta = 0

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
            print('Hubo convergencia, n_iter = ' + str(i+1))
            return x, delta, i+1

    #Si todavia no salio es que no hubo convergencia:
    raise ValueError('No hubo convergencia')
    return x. delta, i+1

#Intervalo para buscar raiz
a = 0.0
b = 2.0

#Parametros para el algoritmo
a_tol1 = 0.5e-15
a_tol2 = 1e-13
n_max = 100



print('----------------')
print('Metodo Newton Raphson')
print('----------------')
print('')
print('Funcion f1, a_tol = '+str(a_tol1))
iterator=0
root=10
data1 = []
print('K,Yk,delta Y,deltaY/Yk,lambda,p')
root, iterator = NewtonRaphson(f1,df1,root,a_tol1,iterator,n_max,data1)
print('raiz = ' +str(root))
print('n_ite= ' +str(iterator))
print('')
# Graficar(root,iterator,trace0descr,trace1descr,graphtitle)
#Graficar(root,iterator,'funcion analizada y = xˆ2 - 2','Funcion y=x','Gráfica de análisis')

print('----------------')
print('Metodo punto fijo')
print('----------------')
print('')
print('Funcion gf1, a_tol = '+str(a_tol1))
iterator=0
root=10
data2 = []
print('K,Yk,delta Y,deltaY/Yk,lambda,p')
root, iterator = FixedPoint(gf1,root,a_tol1,iterator,n_max,data2)
print('raiz = ' +str(root))
print('n_ite= ' +str(iterator))
print('')


print('----------------')
print('Metodo biseccion')
print('----------------')
print('')
print('Funcion f1, a_tol = '+str(a_tol1))
data3 = []
r, delta, n_iter = bisec(f1, a, b, a_tol1, n_max,data3)
print('raiz = ' +str(r))
print('delta= ' +str(delta))
print('n_ite= ' +str(n_iter))
print('')

GraficarDiferencia(data1,data2,data3)


"""
print('----------------')
print('Metodo secante')
print('----------------')
print('')
print('Funcion f1, a_tol = '+str(a_tol1))
r, delta, n_iter = secante(f1, a, b, a_tol1, n_max)
print('raiz = ' +str(r))
print('delta= ' +str(delta))
print('n_ite= ' +str(n_iter))
print('')
print('Funcion f1, a_tol = '+str(a_tol2))
r, delta, n_iter = secante(f1, a, b, a_tol2, n_max)
print('raiz = ' +str(r))
print('delta= ' +str(delta))
print('n_ite= ' +str(n_iter))
print('')

print('----------------')
print('Metodo brent')
print('----------------')
print('')
print('Funcion f1, a_tol por defecto para la funcion')
#https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.brentq.html
r, results = brentq(f1, a, b, full_output=True)
print('raiz = ' +str(r))
print('Resultados: ')
print(results)

"""


