import pandas as pd
import numpy 
import time
print("Modelo UNIFAC by Dawffydd Riv")

#########################	Entrada de variables  

####	<<<<<<<<<<<<<<<<<<<Entrada de parámetros para UNIFAC>>>>>>>>>>>>>>>>>>>>>>>>>>>

compuestos=pd.read_csv('BBDDCompuestos.csv') 

#### Parámetros de las sustancias en la mezcla(en numeros naturales)
print(compuestos.iloc[1:,0:1])
ielegido=int(input("Elige hasta que número será el rango del calculo de la mezcla de compuestos""\n"))  
indice=ielegido+1
df0=compuestos.iloc[1:indice,1:109]#Aqui el 5 puede ser introducido por otro valor creando una variable mediante input 
df0x=df0.values 
matrix0=numpy.array(df0x) 
parametros=matrix0.astype(numpy.float) 
#print(parametros) 
#print("Forma de la matriz parametros""\n",parametros.shape) 
#### Fin de la matriz con la cantidad de coeficientes en numeros naturales

arevol=pd.read_csv('Area y Volumen.csv') 

#### Obtención de la mtriz ri (Volumen molar relativo=vmr)
df1=arevol.iloc[0:108,2:3] 
df1x=df1.values 
matrix1=numpy.array(df1x) 
Ri_1=matrix1.astype(numpy.float) 
#print(Ri_1) 
#print("Forma de la matriz Ri""\n",Ri_1.shape)
ri=numpy.mat(parametros) * numpy.mat(Ri_1) 
#print("Volumenes moleculares Relativos","\n", ri) 
#### Fin de la matriz ri 

#### Obtención de la mtriz qi (Área molar relativa)
df2=arevol.iloc[0:108,3:4] 
df2x=df2.values 
matrix1=numpy.array(df2x) 
Qi_1=matrix1.astype(numpy.float) 
qi=numpy.mat(parametros) * numpy.mat(Qi_1) 
#print("Áreas moleculares Relativas","\n", qi)
#### Fin de la matriz qi (Área molar relativa)

#### Obtención de la matriz Gk=Vk*Qi
qi2=numpy.reshape(Qi_1,(1,108))
Gk=qi2*parametros
#print("La matriz Gk""\n",Gk)
#print("Forma de la matriz Gk""\n",Gk.shape)
#### Fin obtendcion de la matriz Gk

#### Definiendo la temperatura de la mezcla	
Temp=float(input("Introduce la temperatura de la mezcla en grados Kelvin ""\n"))
#### Fin de Temperatura
 #### Calculo energia funcional de la mezcla

#T=334.85
interac=pd.read_csv('interacciones.csv')
dfi=interac.iloc[1:109,2:110]
dfix=dfi.values 
it=numpy.array(dfix) 
itr=it.astype(numpy.float) 
itrx=itr*(-1)
it=itrx/Temp
thau=numpy.exp(it)
#print("Forma de la matriz Thau""\n",thau.shape)
#*print(thau)
#### Fin cálculo de energia funcional de la mezcla

#### Energía "parcial molar" por especie    (4x119)
sk=numpy.mat(Gk) * numpy.mat(thau)
#*print(sk)
#print("Forma de la matriz Sk""\n",sk.shape)
#### Fin del cálculo de la energía "parcial molar" por especie

####	<<<<<<<<<<<<<<<<<<< Fin entrada de constantes para UNIFAC>>>>>>>>>>>>>>>>>>>>>>>>>>>

#### <<<<<<<<<<<<<<Obtención de las solubilidades de referencia >>>>>>>>>>>>>

#### Entalpias de fusion para cálculo de solubilidad
dfe=compuestos.iloc[1:indice,120:121]#Aqui "1:indice" es la entrada de la variable índice
dfent=dfe.values 
matrixent=numpy.array(dfent) 
entalpias=matrixent.astype(numpy.float) 
#print(entalpias)
#### Fin entalpias

#### Temperaturas de fusión para cálculo de solubilidad
dfe=compuestos.iloc[1:indice,121:122]#Aqui "1:indice" es la entrada de la variable índice
dfTx=dfe.values 
matrixTx=numpy.array(dfTx) 
Tempfusion=matrixTx.astype(numpy.float) 
#print(Tempfusion)
#### Fin temperaturas de fusión

ener=(entalpias/(0.008314472*Tempfusion))*((Tempfusion/Temp)-1)
solub_1=numpy.exp(ener)
print("Las solubilidades de referencia para cada componente son""\n",solub_1)
#### <<<<<<<<<<<<<<Fin solubilidades de referencia >>>>>>>>>>>>>

############################### Fin de entrada de variables >>>>>>>>>>>>>>>>>


#### Captura de las concentraciones de cada compuesto de la mezcla
matrixC=[]
print("Introduce las concentraciones de los",ielegido,"compuestos en la mezcla en su orden de captura")
for i in range(ielegido):
		#print("Introduce las concentraciones del compuesto")
		datos=float(input("Da la  concentración  del compuesto con el indice " + str(i+1)+"  "))
		validacion=input("Son correctos los datos?  ")
	
		while validacion =="no":
			print("Vuelva a introducir los datos")
			datos=float(input("Da su cantidad "))
			validacion=input("Son correctos los datos?  ")
			if validacion != "no":
				#matrixC.append(datos)
				break;
		
		matrixC.append(datos)
	
print("\n")	
concentraciones=numpy.array(matrixC)
C0=concentraciones.astype(numpy.float)
print("Las valores de las concentraciones usadas en forma de matriz serán : ""\n", C0)
print("La forma de la matriz de concentraciones es == ",C0.shape)
#### Fin de captura de las concentraciones de cada compuesto de la mezcla

############ Fin de la entrada de variables

########################################	Aquí empezaría el loop	#####################
#>>>>>>>>>>>>>>>
#dif=numpy.zeros([])

#while numpy.all(C0>=0):
####>>>>Empiezan calculos
#### Calculando Ji e Li

##Calculando Ji
c0xri=numpy.mat(C0) * numpy.mat(ri)
#print("Producto de Concentración por volumen molecular relativo" "\n",c0xri)
#print("La forma de la matriz c0xri  ""\n",c0xri)
Ji=numpy.mat(ri)/numpy.mat(c0xri)
#print("Los valores de ji son""\n",Ji)

##Calculando Li
c0_xqi=numpy.mat(C0) * numpy.mat(qi)
#print("Producto de Concentración por area molecular relativa" "\n",c0_xqi)
Li=numpy.mat(qi)/numpy.mat(c0_xqi)
#print("Los valores de ji son""\n",Li)
#### Fin de calculando Ji e Li

#### Teta o area relativa por especie(Gk*xi)
teta=numpy.mat(C0) * numpy.mat(Gk)
#print("La matriz teta""\n",teta)
#print("Forma de la matriz teta""\n",teta.shape)
#### Fin de la matriz Teta

#### Energía de la mezcla por especie   (1x119)
nuk=numpy.mat(C0) * numpy.mat(sk)
#print("Forma de la matriz nuk""\n",nuk.shape)
#### Fin del cálculo de la mezcla por especie

#### Comenzando la parte final LnRi
a1=numpy.log(sk/nuk)
a2=numpy.multiply(Gk,a1)
a22=numpy.array(a2)
a3=numpy.sum(a22,axis=1, keepdims=True)

b1=(sk/nuk)
b2=numpy.multiply(teta,b1)
b22=numpy.array(b2)
b3=numpy.sum(b22,axis=1, keepdims=True)


c1=numpy.log(Li)
c2=1-c1
c21=numpy.multiply(qi,c2)
c22=numpy.array(c21)
c3=numpy.sum(c22,axis=1, keepdims=True)
LnRi=c3-(b3-a3)
#i print("El vector resultante LnRi""\n" , LnRi)
#i print("La forma de esta matriz LnRi ""\n",LnRi.shape)
#### Fin del cálculo de LnRi

#### Comenzando la parte final LnCi
d1=numpy.log(Ji/Li)
d2=Ji/Li
d3=1-d2+d1
d4=5*numpy.multiply(qi,d3)
d5=numpy.log(Ji)
LnCi=1-Ji+d5-d4
#print("El vector resultante LnCi""\n" , LnCi)
#### Fin del cálculo de LnCi

#### Comienzo de Suma de LnRi + LnCi
e=LnRi+LnCi
print("La suma de la parte residual y combinatoria es ""\n",e) 
#### Fin de la suma de LnRi + LnCi

#### Obtención del coeficiente de actividad 
gama=numpy.exp(e)
print("Los valores del coeficiente de actividad de la mezcla son""\n",gama)
#### Fin del cálculo del coeficiente de actividad


##Iniciando codigo prueba para el equilibrio solido-liquido

#### Producto gama por concentraciones
C0_1=numpy.reshape(C0,(ielegido,1))
gamaA=numpy.multiply(C0_1,gama)
print("Las solubilidades de referencia son""\n",solub_1)
print("Impresion de la multiplicacion C0*gama""\n",gamaA)
gaton=numpy.array(gamaA).flatten('F')
#print("Los valores del producto C0 * gama  son""\n",gaton)
#print("La forma de la matriz del producto C0*gama >>>  ",gaton.shape)
	


err=((gamaA-solub_1)/solub_1)
#erro=numpy.absolute(gamaA-solub_1)/solub_1
print("El error calculado para cada componente es""\n",err)
#if numpy.all(err<=0.1):
#	gama_R=1/gama
#	ener_2=(entalpias/(0.008314472*Tempfusion))*((Tempfusion/Temp)-1)
#	C0_resultante=numpy.exp(ener_2)
#	print("Finalmente, la fraccióm mol de los compuestos es""\n",C0_resultante)
#	break;
	#time.sleep(3)
#else:	
#	C0=C0-0.00005
	#time.sleep(3)
