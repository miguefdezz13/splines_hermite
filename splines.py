"""Splines cúbicos paramétricos con factorización LU, diferencias divididas de Hermite y Horner.

Funciones extraídas del notebook splines_hermite.ipynb.
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

def Matrices_LU(A): 
    
    n = len(A) # Tamaño de la matriz A
    
    # Se crean matrices L y U de ceros con el mismo tamaño que A
    L = np.zeros((n, n))
    U = np.zeros((n, n))
    
    # Algoritmo de la variante Doolittle
    for i in range(n):
        L[i, i] = 1
        for j in range(i, n):
                s = sum(L[i, k] * U[k, j] for k in range(i))
                U[i, j] = A[i, j] - s
            
        for j in range(i+1, n):
                s = sum(L[j, k] * U[k, i] for k in range(i))
                L[j, i] = (A[j, i] - s) / U[i, i]
            
    return (L,U)

def Sust_Hacia_Delante(L,Vector): 
    # Sustitución hacia delante Ly = b
    n = len(L) # Tamaño de la matriz L
    y = np.zeros(n) 
    for i in range(n):
            y[i] = Vector[i] - sum(L[i, j] * y[j] for j in range(i))
    return y

def Sust_Hacia_Atrás(U,y):
# Sustitución hacia atrás Ux = y
    n = len(U) # Tamaño de la matriz U
    x = np.zeros(n)
    for i in range(n-1, -1, -1):
        x[i] = (y[i] - sum(U[i, j] * x[j] for j in range(i+1, n))) / U[i, i]
    
    return x

def Resolución_Sistema_LU(A, Vector):
    resultado1 = Matrices_LU(A)
    # Utiliza el resultado de la función Matrices_LU
    y = Sust_Hacia_Delante(resultado1[0], Vector)
    sol = Sust_Hacia_Atrás(resultado1[1], y)
    return sol

def calcula_lista_Deriv(lista_x, lista_f):
    
    # Los puntos son x_0, ..., x_n
    n=len(lista_x)-1
    
    # Calculo los valores de h y de delta
    h = np.zeros(n)
    delta = np.zeros(n)
    
    for i in range(0, n):
        h[i]=lista_x[i+1]-lista_x[i]
        delta[i] = (lista_f[i+1]-lista_f[i]) / h[i]
        
    # Obtengo la matriz de coeficientes y de los respectivos términos independientes
    A = np.zeros((n+1,n+1))
    B = np.zeros(n+1)
    
    # Obtengo en primer lugar la fila 0 y la fila n (son casos distintos)
    # Fila 0
    A[0][0] = 2/h[0]
    A[0][1] = 1/h[0]
    B[0] = delta[0]/h[0]

    # Fila n (ultima)
    A[n][n-1] = 1/h[n-1]
    A[n][n] = 2/h[n-1]
    B[n] = delta[n-1]/h[n-1]
    
    # Obtengo el resto de filas, iterando sobre cada fila
    for i in range(1, n):
        A[i][i-1] = 1/h[i-1]
        A[i][i+1] = 1/h[i]
        A[i][i]   = 2*(A[i][i-1] + A[i][i+1])
        
        B[i] = delta[i]/h[i] + delta[i-1]/h[i-1]
    
    B=3*B
    
    # Resuelvo el sistema mediante un algoritmo basado factorizacion LU y sustitucion hacia atrás
    # y hacia delante
    lista_d = Resolución_Sistema_LU(A, B)
    
    return(lista_d)

def calc_dif_div_Hermite (lista_x, lista_f, lista_d):
    n = len(lista_x) # Numero de valores de las abcisas
    
    dif_divididas = np.zeros([n,n])

    # La tabla se va completando por columnas
    # En primer lugar, completo la columna 0
    for i in range(n):
        dif_divididas[i][0]=lista_f[i]
    
    # Para las demas utilizo los siguientes bucles para rellenarlas
    for k in range(1,n):
        for i in range(0, n-k):
            if (lista_x[i+k] == lista_x[i]):
                # Como se pasan las abcisas duplicadas, se divide entre 2
                dif_divididas[i][k] = lista_d[int(i/2)]
            
            else:
                dif_divididas[i][k] = (dif_divididas[i+1][k-1] - dif_divididas[i][k-1]) / (lista_x[i+k] -lista_x[i])

    # Introducimos los valores de las abcisas a la izquierda de la matriz
    # De esta manera, podemos reproducir exactamente el metodo que utilizabamos en papel
    adicional = np.array(lista_x).reshape(len(lista_x), 1)
    for i in range(n):
        adicional[i]=lista_x[i]
        
    dif_divididas = np.hstack((adicional, dif_divididas))
    return dif_divididas

# Calculador Lista Polinomios
def spline(lista_x, lista_y):
    lista_d = calcula_lista_Deriv(lista_x, lista_y)
    n = len(lista_x) - 1  # Ajustar la longitud de la lista
    Lista_Polinomios = []

    for i in range(n):
        x_inter_x = np.array([lista_x[i], lista_x[i], lista_x[i + 1], lista_x[i + 1]])  # Corregir nombres de variables
        f_inter_x = np.array([lista_y[i], lista_y[i], lista_y[i + 1], lista_y[i + 1]])  # Corregir nombres de variables
        d_inter_x = np.array([lista_d[i], lista_d[i + 1]])  # Corregir nombres de variables
        # Aplico el metodo hermite para calcular el polinomio que se interpola en dicho intervalo
        Lista_Polinomios.append(calc_dif_div_Hermite(x_inter_x, f_inter_x, d_inter_x))  # Añadir polinomio a la lista

    return Lista_Polinomios

def evaluar_Polinomio_Horner_Hermite(tabla_dif_div_Hermite, x_evaluar):
    n = len(tabla_dif_div_Hermite) - 1  # Ajustar el indice para evitar el error de indice fuera de rango
    resultado = tabla_dif_div_Hermite[0][n+1]
    
    # Aplicar método de Horner
    for i in range(n - 1, -1, -1):  # Iterar desde el penultimo coeficiente hacia atras
        resultado = resultado * (x_evaluar - tabla_dif_div_Hermite[i][0]) + tabla_dif_div_Hermite[0][i+1]
        
    return resultado

def Evalua_Intervalo_en_Polinomio(pol_x, array_abcisas, lista_x):
    n = len(array_abcisas)
    array_imagen = np.zeros(n)

    for i in range(n):
        # Con ese -1 me refiero al ultimo de todos
        if(array_abcisas[i]<=lista_x[-1] and array_abcisas[i] >= lista_x[0]):
            # Buscar intervalo en el que se encuentra array_abcisas[i]
            j = 0
       
            while j < len(lista_x) and j < len(pol_x) and array_abcisas[i] >= lista_x[j] :
                j += 1

            # Evaluar polinomio en el punto correspondiente al intervalo
            array_imagen[i] = evaluar_Polinomio_Horner_Hermite(pol_x[j-1], array_abcisas[i])
 
        else: 
            print("Algunos de dichos puntos se salen del intervalo en el que esta definido el spline")

    return array_imagen

def Divide_Matrix_t(lista_x, lista_f):
    n = len(lista_x)
    matrizX = np.empty((n, 2))
    matrizY = np.empty((n, 2))

    for i in range(0,n):
        matrizX[i][1] = lista_x[i]
        matrizX[i][0] = i+1  # Coordenadas de x
        matrizY[i][0] = i+1  # Coordenadas de y
        matrizY[i][1] = lista_f[i]

    return matrizX, matrizY

def plotspline(matriz_puntos, color):
    lista_x = matriz_puntos[:,0]
    lista_f = matriz_puntos[:,1]
    # Divide en dos matrix_puntos
    Array_Matrix_Separado = Divide_Matrix_t(lista_x, lista_f)

    # Hallo hermite y vectores de polinomios
    array_polinomios_X = spline(Array_Matrix_Separado[0][:,0], Array_Matrix_Separado[0][:,1])
    array_polinomios_Y = spline(Array_Matrix_Separado[1][:,0], Array_Matrix_Separado[1][:,1])

    # Inventar array x
    n = len(Array_Matrix_Separado[0])
    abcisas_X = np.linspace(Array_Matrix_Separado[0][0,0], Array_Matrix_Separado[0][n-1,0], 400)
    abcisas_Y = np.linspace(Array_Matrix_Separado[1][0][0], Array_Matrix_Separado[1][n-1][0], 400)

    # Hacer imagenes
    imagenes_X = Evalua_Intervalo_en_Polinomio(array_polinomios_X, abcisas_X,Array_Matrix_Separado[0][:,0])
    imagenes_Y = Evalua_Intervalo_en_Polinomio(array_polinomios_Y,abcisas_Y, Array_Matrix_Separado[1][:,0])

    # Dibujar en el plano
    plt.plot(imagenes_X, imagenes_Y,color)

def Dibujar_Completo(array_matriz_puntos, colores, color_fondo='k'):
    # Crear figura y ejes
    fig, ax = plt.subplots()
    
    # Establecer color de fondo
    ax.set_facecolor(color_fondo)
    
    # Dibujar cada conjunto de puntos
    for i in range(len(array_matriz_puntos)):
        plotspline(array_matriz_puntos[i], colores[i])
        
    ax.set_title('MÉTODOS NUMÉRICOS')
    # Mostrar la gráfica
    plt.show()


