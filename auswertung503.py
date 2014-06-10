'''
Created on 16.04.2014

@author: julian
'''
from numpy import *
from uncertainties import ufloat
from uncertainties.umath import *
import matplotlib.pyplot as plt


def make_LaTeX_table(data,header, flip= 'false', onedim = 'false'):
    output = '\\begin{table}\n\\centering\n\\begin{tabular}{'
    #Get dimensions
    if(onedim == 'true'):
        if(flip == 'false'):
        
            data = array([[i] for i in data])
        
        else:
            data = array([data])
    
    row_cnt, col_cnt = data.shape
    header_cnt = len(header)
    
    if(header_cnt == col_cnt and flip== 'false'):
        #Make Format
        
        for i in range(col_cnt):
            output += 'S'
        output += '}\n\\toprule\n{'+ header[0]
        for i in range (1,col_cnt):
            output += '} &{ ' + header[i]
        output += ' }\\\\\n\\midrule\n'
        for i in data:
            if(isinstance(i[0],(int,float,int32))):
                output += str( i[0] )
            else:
                output += ' ${:L}$ '.format(i[0])
            for j in range(1,col_cnt):
                if(isinstance(i[j],(int,float,int32))):
                    output += ' & ' + str( i[j])
                else:
                    output += ' & ' + str( i[j]).replace('/','')
                
            output += '\\\\\n'
        output += '\\bottomrule\n\\end{tabular}\n\\label{}\n\\caption{}\n\\end{table}\n'
                            
        return output

    else:
        return 'ERROR'



    
def err(data):
    mean = data.mean()
    N = len(data)
    err = 0
    for i in data:
        err += (i - mean)**2
    err = sqrt(err/((N-1)*N))
    return ufloat(mean,err)


def lin_reg(x,y):
    N = len(x)
    sumx = x.sum()
    sumy = y.sum()
    sumxx = (x*x).sum()
    sumxy = (x*y).sum()
    m = (sumxy - sumx*sumy/N)/(sumxx- sumx**2/N)
    b = sumy/N - m*sumx/N
    
    sy = sqrt(((y - m*x - b)**2).sum()/(N-1))
    m_err = sy *sqrt(N/(N*sumxx - sumx**2))
    b_err= m_err * sqrt(sumxx/N)
    return ufloat(m,m_err), ufloat(b,b_err)
    
    
# Auswertung Teil A

# Messwerte einlesen
tnull, tauf, tab, R ,temp= loadtxt("dataA.txt", unpack=True)
s=0.0005

# Erste Tabelle erstellen (Messwerttabelle)
print "Tabelle 1:"
data=array([tnull,tauf,tab,R]).T
print make_LaTeX_table(data,[r'$t_{Null}$[s]', r'$t_{auf}$[s]', r'$t_{ab}$[s]', r'R[M$\omega$]'], flip= 'false', onedim = 'false')
vnull=s/tnull
vauf=s/tauf
vab=s/tab
test1=2*vnull
test2=vab-vauf
test=test2/test1
print test

# gute Messwerte einlesen (Aus zweiter Datei)
tnull, tauf, tab, R , temp, visko= loadtxt("dataA2.txt", unpack=True)
vnull=s/tnull
vauf=s/tauf
vab=s/tab
test1=2*vnull
test2=vab-vauf
test=test2/test1
# Viskositaet in SI umrechnen:
visk=visko*10**(-5)

# Konstanten:
g=9.81
rhoLuft=1.1839
rhoOel=886


# Viskositaet und Radius berechnen:

r = ((9*visk*(vab-vauf))/(2*g*(rhoOel-rhoLuft)))**0.5
print "Radien:"
print r

# Korrektur der Viskose, dabei in SI umgerechnet...
b=(6.17*133.322)*10**(-5)
viskkorr=visk*(1/(1+(b/(101325*r))))

print "Visk und Viskkorrigiert:"
print visk
print viskkorr

# Berechnung der Ladung

q=0.03*pi*viskkorr*((9*viskkorr*(vab-vauf))/(4*g*(rhoOel-rhoLuft)))**(0.5) *((vab+vauf)/(297))
print "Ladungen:"
print q

print "Tabelle 2:"
data=array([vauf,vab,temp,visk,viskkorr,r,q]).T
print make_LaTeX_table(data,[r'$v_{auf}$[m/s]',r'$v_{ab}$[m/s]', r'T[$^\circ C$]', r'$\eta_{L}[10^{-5}Nsm^{-2}]$',r'$\eta[10^{-5}Nsm^{-2}]$',r'r[m]',r'q[C]'], flip= 'false', onedim = 'false')

# Erster Plot
plt.ylim(0,100)
plt.xlabel('gemessene Ladungen [C]')
plt.xlim(0,6.1*10**(-19))
plt.plot(q,[1,1,1,1,1,1,1,1,1,1,1],'.')
plt.show()
plt.savefig("plot1.png")
plt.close()

print "Tabelle 2:"
data=array([vauf,vab,temp,visk,viskkorr,r,q]).T
print make_LaTeX_table(data,[r'$v_{auf}$[m/s]',r'$v_{ab}$[m/s]', r'T[$^\circ C$]', r'$\eta_{L}[10^{-5}Nsm^{-2}]$',r'$\eta[10^{-5}Nsm^{-2}]$',r'r[m]',r'q[C]'], flip= 'false', onedim = 'false')

# Auswertung Teil B

# Messwerte einlesen
tab, tauf, R, temp, visko= loadtxt("dataC.txt", unpack=True)
s=0.0005

# Zweite Tabelle erstellen (Messwerttabelle)
print "Tabelle 3:"
data=array([tauf,tab,R]).T
print make_LaTeX_table(data,[r'N', r't[s]', r'$\sigma_{I,rel}$[\%]'], flip= 'false', onedim = 'false')

vauf=s/tauf
vab=s/tab

# Viskositaet in SI umrechnen:
visk=visko*10**(-5)



# Viskositaet und Radius berechnen:

r = ((9*visk*(vab-vauf))/(2*g*(rhoOel-rhoLuft)))**0.5
print "Radien:"
print r

# Korrektur der Viskose, dabei in SI umgerechnet
b=(6.17*133.322)*10**(-5)
viskkorr=visk*(1/(1+(b/(101325*r))))



# Berechnung der Ladung

q=0.03*pi*viskkorr*((9*viskkorr*(vab-vauf))/(4*g*(rhoOel-rhoLuft)))**(0.5) *((vab+vauf)/(150))
print "Ladungen:"
print q

# Zweiter Plot
plt.ylim(0,100)
plt.xlabel('gemessene Ladungen [C]')
plt.xlim(0,9*10**(-19))
plt.plot(q,[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],'.')
plt.show()
plt.savefig("plot2.png")
plt.close()

print "Tabelle 4:"
data=array([vauf,vab,temp,visk,viskkorr,r,q]).T
print make_LaTeX_table(data,[r'$v_{auf}$[m/s]',r'$v_{ab}$[m/s]', r'T[$^\circ C$]', r'$\eta_{L}[10^{-5}Nsm^{-2}]$',r'$\eta[10^{-5}Nsm^{-2}]$',r'r[m]',r'q[C]'], flip= 'false', onedim = 'false')




print "AUSWERTUNG TEIL C:"
#Auswertung Teil C
# Messwerte einlesen
tnull, tauf, tab, R, temp, visko= loadtxt("dataB.txt", unpack=True)
s=0.0005

Tab=(tab[0]+tab[1]+tab[2])/3
Tauf=(tauf[0]+tauf[1]+tauf[2])/3
vab=s/Tab
vauf=s/Tauf

visk=visko[0]*10**(-5)

r = ((9*visk*(vab-vauf))/(2*g*(rhoOel-rhoLuft)))**0.5

# Korrektur der Viskose, dabei in SI umgerechnet...
b=(6.17*133.322)*10**(-5)
viskkorr=visk*(1/(1+(b/(101325*r))))

q=0.03*pi*viskkorr*((9*viskkorr*(vab-vauf))/(4*g*(rhoOel-rhoLuft)))**(0.5) *((vab+vauf)/(150))
print "Ladungen:"
print q

'''
# Erster Plot
plt.xlim(290,710)
plt.ylim(-1,62)
plt.xlabel('Zaehlrohrspannung U[V]')
plt.ylabel('Impulsrate I [1/s]')
plt.errorbar(u,I,yerr=sigmaI,xerr=None,fmt=None,ecolor='red')
plt.plot(u,I,'.')
plt.show()
plt.savefig("plot1.png")
plt.close()

# Zweiter Plot (Regressionsgerade)

# Arrays mit benoetigten Werten erstellen:
u2=u[2:] 
I2=I[2:]
sigmaI2=sigmaI[2:]

# Lineare Regression
m,b = lin_reg(u2,I2)

print "M:"
print m
print "B:"
print b

plt.xlim(335,710)
plt.xlabel('Zaehlrohrspannung U[V]')
plt.ylabel('Impulsrate I [1/s]')
plt.errorbar(u2,I2,yerr=sigmaI2,xerr=None,fmt=None,ecolor='red')
plt.plot(u2,I2,'.')
plt.plot(u2,u2*m.n+b.n)
plt.show()
plt.savefig("plot2.png")
plt.close()
'''


