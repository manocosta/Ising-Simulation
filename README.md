# Ising-Simulation
/*
		Este programa se realizo con el objetivo de poder realizar simulaciones del modelo de 
	Ising mediante el metodo de Metropolis Montecarlo permitiendo definir los parametros de
	estas simulaciones en el llamado al programa. 
	
	Un posible llamado al programa es el siguiente:
				./programa.exe SEED 260572 L 32 INIT_T 5 B 0.5
		En este caso se setea la semilla (SEED = 260572) que se utilizara para iniciar el proceso c
	omo, el tamaño de la red que se utilizara (L = 32 es decir una red de 32x32), la temperatura 
	inicial (INIT_T = 5) y el campo magnetico (B = 0.5). Como se puede ver, cada parametro se setea 
	a traves de escribir el comando correspondiente al parametro que se desea setear seguido del 
	valor deseado.
		A continuacion podemos ver los parametros que se pueden setear y los comandos 
	correspondientes a cada uno de estos:
		  SEED   : Semilla utilizada para iniciar el programa.
		    L    : Tamaño de la red (LxL).
		 INIT_T  : Temperatura inicial del sistema (dado en relacion a la k de Boltzmann).
		 STOP_T  : Temperatura en que se finalizara el ciclo.
		INIT_HOT : 1 si se inicializa la red con un sistema caliente (random spins) o 0 si se inicializa frio (todos los espines iguales).
		   D_T   : Tamaño de paso al cambiar de temperatura.
		    J    : Constante de la interaccion de Ising entre los espines vecinos.
			B    : Campo magnetico externo.
		  PROB   : Probabilidad con la que setearan los spines en caso de iniciar el sistema caliente.
		K_DESCOR : Cantidad de pasos utilizados para descorrelacionar mediciones.
		K_TERMAL : Cantidad de pasos requeridos para termalizar el sistema luego de modificar la temperatura.
		 N_MEDS  : Numero de mediciones que se promediaran para cada temperatura.
		  PATH   : Directorio en que se desean guardar los resultados.
		  FILE   : Nombre del archivo en que se guardaran los resultados.
*/
