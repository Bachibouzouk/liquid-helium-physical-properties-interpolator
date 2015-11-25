# -*- coding: utf-8 -*-
"""
@author: Pierre-Francois Duc pfduc@physics.mcgill.ca

@author: pfduc
"""
import numpy as np
import pylab as plt
from scipy import interpolate
from os.path import join

Mm=4.002602e-3 #kg/mol from http://webbook.nist.gov/cgi/cbook.cgi?ID=C7440597&Type=JANAFG&Plot=on
Na=6.02214129*1e23 #1/mol 
R=8.3144621 #J/(mol K) (±0.0000075) from http://physics.nist.gov/cgi-bin/cuu/Value?r
Rs=R/Mm  #J/(kg K)

kb=1.3806488*1e-23 #J K-1
h_planck=6.62606957*1e-34 #J s
R=8.3144621 #J/(mol K) (±0.0000075) from http://physics.nist.gov/cgi-bin/cuu/Value?r
Rs=R/Mm  #J/(kg K)
m_He4=Mm/Na


#This is the path to the data sets
he4_data_path="expt_datasets/"

def density(P,T):
    """
    Pressure in psi and temperature in K, this fonction interpolates between many different sources of helium properties (mainly Donnelly)
    Depending on the pressure or temperature range, this could return results with more or less errors.
    The most sensitive approximations are around the superfluid transition temperature T_lambda    
    """
    answer=np.array([])
    if np.size(P)>1 and np.size(T)>1:
        answer=[]
        for t in T:
            answer_line=np.array([])
            for p in P:
                answer_line=np.append(answer_line,density(p,t))
            answer.append(answer_line)
        return np.vstack(answer)
    elif np.size(T)>1 and np.size(P)==1:
        for t in T:
            answer=np.append(answer,density(P,t))
        return answer
    elif np.size(P)>1 and np.size(T)==1:
        for p in P:
            answer=np.append(answer,density(p,T))
        return answer
    elif np.size(P)==1 and np.size(P)==np.size(T):
        return np.squeeze(interpolate_density(P,T))
    else:
        print "error in He4Property.density, T and P arrays are empty"
    

def normal_density(P,T):
    """
    Pressure in psi and temperature in K, this fonction interpolates between many different sources of helium properties (mainly Donnelly)
    Depending on the pressure or temperature range, this could return results with more or less errors.
    The most sensitive approximations are around the superfluid transition temperature T_lambda    
    """
    answer=np.array([])
    if np.size(P)>1 and np.size(T)>1:
        answer=[]
        for t in T:
            answer_line=np.array([])
            for p in P:
                answer_line=np.append(answer_line,normal_density(p,t))
            answer.append(answer_line)
        return np.vstack(answer)
    elif np.size(T)>1 and np.size(P)==1:
        for t in T:
            answer=np.append(answer,normal_density(P,t))
        return answer
    elif np.size(P)>1 and np.size(T)==1:
        for p in P:
            answer=np.append(answer,normal_density(p,T))
        return answer
    elif np.size(P)==1 and np.size(P)==np.size(T):
        if T<=T_lambda(P):
            return (1-superfluid_density_ratio(P,T))*density(P,T)
        else:
            return density(P,T)
    else:
        print "error in He4Property.normal_density, T and P arrays are empty"

    
    
        
def superfluid_density(P,T):
    """
    Pressure in psi and temperature in K, this fonction interpolates between many different sources of helium properties (mainly Donnelly)
    Depending on the pressure or temperature range, this could return results with more or less errors.
    The most sensitive approximations are around the superfluid transition temperature T_lambda    
    """
    answer=np.array([])
    if np.size(P)>1 and np.size(T)>1:
        answer=[]
        for t in T:
            answer_line=np.array([])
            for p in P:
                answer_line=np.append(answer_line,superfluid_density(p,t))
            answer.append(answer_line)
        return np.vstack(answer)
    elif np.size(T)>1 and np.size(P)==1:
        for t in T:
            answer=np.append(answer,superfluid_density(P,t))
        return answer
    elif np.size(P)>1 and np.size(T)==1:
        for p in P:
            answer=np.append(answer,superfluid_density(p,T))
        return answer
    elif np.size(P)==1 and np.size(P)==np.size(T):
        if T<=T_lambda(P):
            return superfluid_density_ratio(P,T)*density(P,T)
        else:
            return 0
    else:
        print "error in He4Property.superfluid_density, T and P arrays are empty"
    

def superfluid_density_ratio(P,T):
    answer=np.array([])
    if np.size(P)>1 and np.size(T)>1:
        answer=[]
        for t in T:
            answer_line=np.array([])
            for p in P:
                answer_line=np.append(answer_line,superfluid_density_ratio(p,t))
            answer.append(answer_line)
        return np.vstack(answer)
    elif np.size(T)>1 and np.size(P)==1:
        for t in T:
            answer=np.append(answer,superfluid_density_ratio(P,t))
        return answer
    elif np.size(P)>1 and np.size(T)==1:
        for p in P:
            answer=np.append(answer,superfluid_density_ratio(p,T))
        return answer
    elif np.size(P)==1 and np.size(P)==np.size(T):
            return np.squeeze(interpolate_density(P,T,fname=join(he4_data_path,"He_superfluid_ratio.dat")))
    else:
        print "error in He4Property.density, T and P arrays are empty"

def interpolate_density(P,T,fname=join(he4_data_path,"He_liquid_density_tables.dat"),test_plot=False):
    if T<0.65:
        print "the data in He_liquid_SVP.dat do not go lower that 0.65K so we cannot get the density for %.4f"%T
    dat=np.loadtxt(fname)
    Tref=dat[:,0]
    Pref=36.7398719 #this is in psi
    dat=np.transpose(dat)
    #selects only the range over which my experiments were conducted, ie less than 2.5at
    dat=dat[1:3]
    dens_fixed_P=[]

    dat_svp=np.loadtxt(fname=join(he4_data_path,"He_liquid_SVP.dat"))

    if T in Tref:
        dens=dat[:,Tref==T]
        SVP=np.squeeze(dat_svp[dat_svp[:,0]==T,1]/6894.7572798677)
        if P<SVP:
            print "He4Properties.interpolate_density : The pressure %.3f is lower than the saturated vapor pressure(%.3f) for the temperature %.3f, it is not liquid anymore"%(P,SVP,T)
        interp1d_over_P=interpolate.interp1d(np.array([SVP,Pref]),np.squeeze(dens),kind='linear')
        if test_plot:
            #    Test plot to see if the interpolation gives a good estimate
            plt.plot(T,interp1d_over_P(P),'-ok')
            for d in dat:
                plt.hold(True)
                plt.plot(Tref,d,'b')
            plt.show()
        try:
            return interp1d_over_P(P)
        except ValueError:
            return np.nan
    else: 
        Tref2=[]
        for t in Tref:
            dens=dat[:,Tref==t]
            SVP=np.squeeze(dat_svp[dat_svp[:,0]==t,1]/6894.7572798677)
            if P>SVP:
                interp1d_over_P=interpolate.interp1d(np.array([SVP,Pref]),np.squeeze(dens),kind='linear')
                try:
                    dens_fixed_P.append(interp1d_over_P(P))
                    Tref2.append(t)
                except ValueError:
                    print "the value of pressure %.2f atm is out of interpolation range (%.2f-%.2f) for the selected temperature %.3f K, please input a value in the following list:"%(P,SVP,Pref,T)
                    print Pref
                
        #interpolate over the interpolated fixed pressure density versus temperature curve    
        interp1d_over_T=interpolate.InterpolatedUnivariateSpline(Tref2,dens_fixed_P)
        if test_plot:
            #    Test plot to see if the interpolation gives a good estimate
            plt.plot(Tref,dens_fixed_P,'-xr')
            for d in dat:
                plt.hold(True)
                plt.plot(Tref,d,'b')
        return interp1d_over_T(T)

def viscosity(P,T):
    answer=np.array([])
    if np.size(P)>1 and np.size(T)>1:
        answer=[]
        for t in T:
            answer_line=np.array([])
            for p in P:
                print p,t
                answer_line=np.append(answer_line,viscosity(p,t))
            answer.append(answer_line)
        return np.vstack(answer)
    elif np.size(T)>1 and np.size(P)==1:
        for t in T:
            answer=np.append(answer,viscosity(P,t))
        return answer
    elif np.size(P)>1 and np.size(T)==1:
        for p in P:
            answer=np.append(answer,viscosity(p,T))
        return answer
    elif np.size(P)==1 and np.size(P)==np.size(T):
        return np.squeeze(interpolate_viscosity(P,T,join(he4_data_path,"He_liquid_viscosity.dat")))
    else:
        print "error in He4Property.viscosity, T and P arrays are empty"
    


#
def interpolate_viscosity(P,T,fname=join(he4_data_path,"He_liquid_viscosity.dat"),test_plot=False):
    """
    this one is first making a linear interpolation over the pressure between all 
    the same temp points, therefore obtaining a curve for a constant pressure, 
    then this curve is itself cubic interpolated to get the temperature dependance.
    """    
    if T<1.2:
        print "the data in He_liquid_viscosity.dat do not go lower that 1.2K so we cannot get the viscosity for %.4f"%T
        return np.nan
    else:
        dat=np.loadtxt(fname)
        Tref=dat[:,0]
        Pref=64.6621746#value is in psi
        dat=np.transpose(dat)
        #selects only the range over which experiments were conducted, ie less than 2.5at
        dat=dat[1:3]
        visc_fixed_P=[]

        dat_svp=np.loadtxt(fname=join(he4_data_path,"He_liquid_SVP.dat"))
    #    print dat_svp[:,0]
        if T in Tref:
            visc=dat[:,Tref==T]
            SVP=np.squeeze(dat_svp[dat_svp[:,0]==T,1]/6894.7572798677)
            if P<SVP:
                print "He4Properties.interpolate_density : The pressure %.3f is lower than the saturated vapor pressure(%.3f) for the temperature %.3f, it is not liquid anymore"%(P,SVP,T)
    
            interp1d_over_P=interpolate.interp1d(np.array([SVP,Pref]),np.squeeze(visc),kind='linear')
            if test_plot:
                #    Test plot to see if the interpolation gives a good estimate
                plt.plot(T,interp1d_over_P(P),'-ok')
                for d in dat:
                    plt.hold(True)
                    plt.plot(Tref,d,'b')
                plt.show()
            try:
                return interp1d_over_P(P)
            except ValueError:
                return np.nan
        else: 

            Tref2=[]
            for t in Tref:
                visc=dat[:,Tref==t]
                SVP=np.squeeze(dat_svp[dat_svp[:,0]==t,1]/6894.7572798677)
                if P>SVP:
                    interp1d_over_P=interpolate.interp1d(np.array([SVP,Pref]),np.squeeze(visc),kind='linear')
                    try:
                        visc_fixed_P.append(interp1d_over_P(P))
                        Tref2.append(t)
                    except ValueError:
                        print "the value of pressure %.2f atm is out of interpolation range (%.2f-%.2f) for the selected temperature %.3f K, please input a value in the following list:"%(P,SVP,Pref,T)
                        print Pref
                    
            #interpolate over the interpolated fixed pressure density versus temperature curve        
            interp1d_over_T=interpolate.InterpolatedUnivariateSpline(Tref2,visc_fixed_P)
            if test_plot:
                #    Test plot to see if the interpolation gives a good estimate
                plt.plot(Tref,visc_fixed_P,'-xr')
                #    plt.plot(T,interp1d_over_T(T),'ok')
                for d in dat:
                    plt.hold(True)
                    plt.plot(Tref,d,'b')
            return interp1d_over_T(T)

def T_lambda(P):
    """
    this function gives the temperature at which the liquid HeI becomes liquid HeII or superfluid for a given pressure
    The argrument can also be iterable.
    """
    if np.size(P)>1:
        answer=np.array([])
        for p in P:
            answer=np.append(answer,T_lambda(p))
        return answer
    else:
        data=np.loadtxt(join(he4_data_path,"Lambda_line.dat"))
        my_interp=interpolate.interp1d(np.squeeze(data[:,0]),np.squeeze(data[:,1]),kind='linear')
        if (data<P).any()==False:
            return data[0,-1]
        else:
            return my_interp(P)
            
def saturated_vapor_pressure(T):
    """
    this function gives the temperature at which the liquid HeI becomes liquid HeII or superfluid for a given pressure
    The argrument can also be iterable.
    The returned value is in psi
    """
    if np.size(T)>1:
        answer=np.array([])
        for t in T:
            answer=np.append(answer,saturated_vapor_pressure(t))
        return answer
    else:
        data=np.loadtxt(fname=join(he4_data_path,"He_liquid_SVP.dat"))
        my_interp=interpolate.interp1d(np.squeeze(data[:,0]),np.squeeze(data[:,1])/6894.7572798677,kind='linear')
        return my_interp(T)
        
def uncertainty(func,P,dP,T,dT):
    """
    this function returns and uncertainty on the value returned by "func" depending on the pressure, temperature
    and your uncertainty in them. As The viscosity and density are changing fast close to the transition temperature
    the actual values of pressure and temperature are as important as the dP and dT.
    """
    answer=np.array([])
    
    #Treating various input senarii
    if np.size(P)>1 and np.size(T)==np.size(P):
        for t,p in zip(T,P):
            answer=np.append(answer,uncertainty(func,p,dP,t,dT))
        return answer
        
    elif np.size(T)>1 and np.size(P)==1:
        for t in T:
            answer=np.append(answer,uncertainty(func,P,dP,t,dT))
        return answer
        
    elif np.size(P)>1 and np.size(T)==1:
        for p in P:
            answer=np.append(answer,uncertainty(func,p,dP,T,dT))
        return answer
        
    elif np.size(P)>1 and np.size(T)>1:
        #make a recursive call to the function
        answer=[]
        for t in T:
            answer_line=np.array([])
            for p in P:
                answer_line=np.append(answer_line,uncertainty(func,p,dP,t,dT))
            answer.append(answer_line)
            
        return np.vstack(answer)
        
    elif np.size(P)==1 and np.size(T)==1:
        if func.__name__=="superfluid_density" and T>T_lambda(P):
            return 0
        else:
            value=func(P,T)
            Prange=[P-dP,P+dP]
            Trange=[T-dT,T+dT]
            #gets all the values of the func when varying P and T within their uncertainties range
            value_uncertainties=[func(p,t) for p in Prange for t in Trange]
            
            low_bound=np.max(value-value_uncertainties)
            up_bound=np.max(value_uncertainties-value)

            return (low_bound+up_bound)/(2*np.sqrt(3))# according to http://www.bipm.org/en/publications/guides/gum the GUM chapter 4.3.8
    else:
        print "error in He4Property.uncertainty, T and P arrays are empty"


if __name__=="__main__":
    
    print viscosity(30,2.3)
    print "The density of liquid helium4 is at P= and T= is equal to %.3f kgm^-3"%(density(30,2.3))
 