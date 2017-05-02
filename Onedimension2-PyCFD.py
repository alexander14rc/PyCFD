class Onedimension:
    """docstring for Onedimension"""
    def __init__(self,nx ,nt,Lx):
        self.nx = nx
        self.nt = nt
        self.Lx = Lx


class Funcion(Onedimension):
    def __init__(self,nx = 41 ,nt = 25,Lx = 2.,dt = 0.025, c = 1):
        super().__init__(nx,nt,Lx)
        self.dt = dt
        self.c = c
        
        
    def wave_lineal(self,nx = 41 ,nt = 25,Lx = 2.,dt = 0.025, c = 1):
        super().__init__(nx,nt,Lx)
        self.dt = dt
        self.c = c
        import matplotlib.pyplot as plt
        from drawnow import drawnow
        import numpy as np
        import time, sys  
        import copy 

  
        def makeFig():
            plt.plot(xList,yList,'b')
            plt.plot(1,Lx)
            plt.xlabel('X (m)')
            plt.ylabel('Y (m)')
            plt.title("ECUACION DE LA ONDA\n Lineal \n Tiempo {:.3f} s".format(inter))
            plt.grid(True)
            plt.annotate(s=u"Celeridad = {} m/s".format(self.c), xy=(0.95, 1.8), xytext=(0.2, 1.8), arrowprops=dict(arrowstyle = "->"))
            plt.savefig("Onda_lineal.png") #Grafico de la Onda en PNG
            plt.savefig("Onda_lineal.jpg") #Grafico de la Onda en JPG
            plt.savefig("Onda_lineal.pdf") #Grafico de la Onda en PDF
            plt.show()
            


        plt.ion() 
        fig=plt.figure() 

       
        dx = Lx/(nx -1)
        u = np.ones(nx)
        A = np.zeros((nt, nx))
        u[int(.5/dx) : int(1./dx + 1)] = 2 # Para valores dentro de un array siempre poner enteros
        # u0 = copy.copy(u)
        un = np.ones(nx)
        for n in range(nt):
            un[:] = u[:]
            for i in range(1,nx):
                u[i] = un[i] - c * (dt/dx)*(un[i]-un[i-1])
            A[n,:] = u
    
        for i in range(nt):
            xList = np.linspace(0,Lx,nx)
            yList = A[i,:]
            inter = i * dt
            drawnow(makeFig)
            plt.pause(0.01)
        drawnow(makeFig)
        plt.show('hold')

    def wave_no_lineal(self,nx = 41 ,nt = 25,Lx = 2.,dt = 0.025):
        super().__init__(nx,nt,Lx)
        self.dt = dt
        import matplotlib.pyplot as plt
        from drawnow import drawnow
        import numpy as np
        import time, sys  
        import copy 

  
        def makeFig():
            plt.plot(xList,yList)
            plt.plot(1,Lx)
            plt.xlabel('X (m)')
            plt.ylabel('Y (m)')
            plt.title("ECUACION DE LA ONDA\n No Lineal \n Tiempo{:.3f} s".format(inter))
            plt.grid(True)
            plt.annotate(s=u"Celeridad ", xy=(0.95, 1.8), xytext=(0.2, 1.8), arrowprops=dict(arrowstyle = "->"))
            plt.savefig("Onda_no_lineal.png") #Grafico de la Onda en PNG
            plt.savefig("Onda_no_lineal.jpg") #Grafico de la Onda en JPG
            plt.savefig("Onda_no_lineal.pdf") #Grafico de la Onda en PDF
            plt.show()


        plt.ion() 
        fig=plt.figure() 

       
        dx = Lx/(nx -1)
        u = np.ones(nx)
        A = np.zeros((nt, nx))
        u[int(.5/dx) : int(1./dx + 1)] = 2 # Para valores dentro de un array siempre poner enteros
        # u0 = copy.copy(u)
        un = np.ones(nx)
        for n in range(nt):
            un[:] = u[:]
            for i in range(1,nx):
                u[i] = un[i] - un[i] * (dt/dx)*(un[i]-un[i-1])
            A[n,:] = u
    
        for i in range(nt):
            xList = np.linspace(0,Lx,nx)
            yList = A[i,:]
            inter = i*dt
            drawnow(makeFig)
            plt.pause(0.1)
        drawnow(makeFig)
        plt.show('hold')


    def heat(self,nx = 41 ,nt = 20,Lx = 2.,sigma = 0.2,nu = 0.3,n0 = 2):
        super().__init__(nx,nt,Lx)
        self.n0 = n0
        self.sigma = sigma
        self.nu = nu
        # nu = valor de la viscosity
        # sigma = es un parametro
        import matplotlib.pyplot as plt
        from drawnow import drawnow
        import numpy as np
        import time, sys  
        import copy 

  
        def makeFig():
            plt.plot(xList,yList,'b')
            plt.plot(1,Lx)
            plt.xlabel('X (m)')
            plt.ylabel('Y (m)')
            plt.title("ECUACION DEL CALOR\n No Lineal \n Tiempo {:.3f} s".format(inter))
            plt.grid(True)
            plt.annotate(s=u"Calor ", xy=(0.95, 1.8), xytext=(0.2, 1.8))
            plt.savefig("Calor_no_lineal.png") #Grafico de la Onda en PNG
            plt.savefig("Calor_no_lineal.jpg") #Grafico de la Onda en JPG
            plt.savefig("Calor_no_lineal.pdf") #Grafico de la Onda en PDF
            plt.show()


        plt.ion() 
        fig=plt.figure() 

       
        dx = Lx/(nx -1)
        u = np.ones(nx)
        dt = sigma*dx**2/nu
        A = np.zeros((nt, nx))
        u[int(0.5/dx) : int(1./dx + 1)] = n0 # Para valores dentro de un array siempre poner enteros
        
        un = np.ones(nx)

        for n in range(nt):
            un[:] = u[:]
            for i in range(1,nx-1):
                u[i] = un[i] + nu * dt/dx**2*(un[i+1]-2*un[i]+un[i-1])
            A[n,:] = u
    
        for i in range(nt):
            xList = np.linspace(0,Lx,nx)
            yList = A[i,:]
            inter = i*dt
            drawnow(makeFig)
            plt.pause(0.1)
        drawnow(makeFig)
        plt.show('hold')


# EJEMPLOS:
#b = Funcion().wave_lineal()
#m = Funcion().wave_lineal(nx = 41 ,nt = 150,Lx = 10,dt = 0.025, c = 1)
#m = Funcion().wave_lineal(41, 150, 10, 0.025,1)
#b = Funcion().wave_no_lineal()
#b = Funcion().wave_no_lineal(nx = 41 ,nt = 150,dt = 0.025,Lx = 10.)
#b = Funcion().wave_no_lineal(41 , 150, 0.025, 10.)
#b = Funcion().heat()
#b = Funcion().heat(nx = 41 , nt = 20, Lx = 2, sigma = 0.2, nu = 0.3, n0 = 2)