import numpy
from dadi import Numerics, PhiManip, Integration
from dadi.Spectrum_mod import Spectrum

'''
Models for testing various three population scenarios.

############################################

Dan Portik
daniel.portik@uta.edu
March 2017
'''

##########################################################################################
#Basic models of (no gene flow / gene flow) between (all / some) population pairs
##########################################################################################
# pops 1,2, and 3 are ECO,SW,SSAN
def model22_split_nomig(params, ns, pts): # OK
    """
    Model with split between pop 1 and (2,3), then split between 2 and 3.
    Migration does not occur between any population pair.
    """
    #6 parameters	
    nu1, nuA, nu2, nu3, T1, T2 = params

    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T1, nu1=nu1, nu2=nuA, m12=0, m21=0)

    phi = PhiManip.phi_2D_to_3D_split_2(xx, phi)
    
    phi = Integration.three_pops(phi, xx, T2, nu1=nu1, nu2=nu2, nu3=nu3, m12=0, m21=0, m23=0, m32=0, m13=0, m31=0)

    fs = Spectrum.from_phi(phi, ns, (xx,xx,xx))
    return fs


def model23_split_symmig_all(params, ns, pts): # OK
    """
    Model with split between pop 1 and (2,3), then split between 2 and 3.
    Migration is symmetrical between all population pairs (ie 1<->2, 2<->3, and 1<->3).
    """
    #10 parameters
    nu1, nuA, nu2, nu3, mA, m1, m2, m3, T1, T2 = params

    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T1, nu1=nu1, nu2=nuA, m12=mA, m21=mA)

    phi = PhiManip.phi_2D_to_3D_split_2(xx, phi)
    
    phi = Integration.three_pops(phi, xx, T2, nu1=nu1, nu2=nu2, nu3=nu3, m12=m1, m21=m1, m23=m2, m32=m2, m13=m3, m31=m3)

    fs = Spectrum.from_phi(phi, ns, (xx,xx,xx))
    return fs

def model24_split_assymmig_all(params, ns, pts): # OK
    """
    Model with split between pop 1 and (2,3), then split between 2 and 3.
    Migration is symmetrical between for-eco, by asym for for-for.
    """
    #10 parameters
    nu1, nuA, nu2, nu3, mA, m1, m2, m3, m4, T1, T2 = params

    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T1, nu1=nu1, nu2=nuA, m12=mA, m21=mA)

    phi = PhiManip.phi_2D_to_3D_split_2(xx, phi)

    phi = Integration.three_pops(phi, xx, T2, nu1=nu1, nu2=nu2, nu3=nu3, m12=m1, m21=m1, m23=m2, m32=m3, m13=m4, m31=m4)

    fs = Spectrum.from_phi(phi, ns, (xx,xx,xx))
    return fs

    
##########################################################################################
#Various models based on ancient migration and contemporary isolation
##########################################################################################

def model25_ancmig_3(params, ns, pts): # OK
    """
    Model with split between pop 1 and (2,3), with gene flow, which then stops. Split 
    between pops 2 and 3, gene flow does not occur at all.
    longest isolation
    """
    #8 parameters
    nu1, nuA, nu2, nu3, mA, T1a, T1b, T2 = params

    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T1a, nu1=nu1, nu2=nuA, m12=mA, m21=mA)
    
    phi = Integration.two_pops(phi, xx, T1b, nu1=nu1, nu2=nuA, m12=0, m21=0)

    phi = PhiManip.phi_2D_to_3D_split_2(xx, phi)
    
    phi = Integration.three_pops(phi, xx, T2, nu1=nu1, nu2=nu2, nu3=nu3, m12=0, m21=0, m23=0, m32=0, m13=0, m31=0)
    
    fs = Spectrum.from_phi(phi, ns, (xx,xx,xx))
    return fs

def model26_ancmig_2(params, ns, pts): # OK
    """
    Model with split between pop 1 and (2,3), with gene flow. Split 
    between pops 2 and 3, and all gene flow ceases.
    shorter isolation
    """
    #7 parameters
    nu1, nuA, nu2, nu3, mA, T1, T2 = params

    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T1, nu1=nu1, nu2=nuA, m12=mA, m21=mA)

    phi = PhiManip.phi_2D_to_3D_split_2(xx, phi)
    
    phi = Integration.three_pops(phi, xx, T2, nu1=nu1, nu2=nu2, nu3=nu3, m12=0, m21=0, m23=0, m32=0, m13=0, m31=0)
    
    fs = Spectrum.from_phi(phi, ns, (xx,xx,xx))
    return fs

def model27_ancmig_1(params, ns, pts): # OK
    """
    Model with split between pop 1 and (2,3), with gene flow. Split 
    between pops 2 and 3 with gene flow, then all gene flow ceases.
    shortest isolation, MODIFIED SO GENE FLOW BETEWEEN 1 and 3
    . gene flow between two forest blocks assymetric, for-eco sym
    """
    #11 parameters
    nu1, nuA, nu2, nu3, mA, m1, m2, m3, m4, T1, T2, T3 = params

    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T1, nu1=nu1, nu2=nuA, m12=mA, m21=mA)
    
    phi = PhiManip.phi_2D_to_3D_split_2(xx, phi)
    
    phi = Integration.three_pops(phi, xx, T2, nu1=nu1, nu2=nu2, nu3=nu3, m12=m1, m21=m1, m23=m2, m32=m3, m13=m4, m31=m4)

    phi = Integration.three_pops(phi, xx, T3, nu1=nu1, nu2=nu2, nu3=nu3, m12=0, m21=0, m23=0, m32=0, m13=0, m31=0)

    fs = Spectrum.from_phi(phi, ns, (xx,xx,xx))
    return fs


def model28_ancmig_4(params, ns, pts): # OK
    """
    Model with split between pop 1 and (2,3), with gene flow. Split
    between pops 2 and 3 with gene flow, with ongoing gene flow
    between 1 and 2/3, then cessation of forest-ecotone gene flow
    ancient for-eco divergence, cessation, but between forest block gene flow
    for-eco symmetric, for-for assymetric
    """
    #11 parameters
    nu1, nuA, nu2, nu3, mA, m1, m2, m3, m4, T1, T2, T3 = params

    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T1, nu1=nu1, nu2=nuA, m12=mA, m21=mA)
    
    phi = PhiManip.phi_2D_to_3D_split_2(xx, phi)
    
    phi = Integration.three_pops(phi, xx, T2, nu1=nu1, nu2=nu2, nu3=nu3, m12=m1, m21=m1, m23=m2, m32=m3, m13=m4, m31=m4)

    phi = Integration.three_pops(phi, xx, T3, nu1=nu1, nu2=nu2, nu3=nu3, m12=0, m21=0, m23=m2, m32=m3, m13=0, m31=0)

    fs = Spectrum.from_phi(phi, ns, (xx,xx,xx))
    return fs


def model29_ancmig_5(params, ns, pts): # OK
    """
    Model with split between pop 1 and (2,3), with gene flow. which
    then ceases prior to 2-3 split. Split between pops 2 and 3 with gene flow, 
    ancient for-eco divergence, cessation pre forest split, but between forest block gene flow
    for-eco sym, for-for asym
    """
    #9 parameters
    nu1, nuA, nu2, nu3, mA, m1, m2, T1a, T1b, T2  = params

    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)


    phi = Integration.two_pops(phi, xx, T1a, nu1=nu1, nu2=nuA, m12=mA, m21=mA)

    phi = Integration.two_pops(phi, xx, T1b, nu1=nu1, nu2=nuA, m12=0, m21=0)

    phi = PhiManip.phi_2D_to_3D_split_2(xx, phi)
    
    phi = Integration.three_pops(phi, xx, T2, nu1=nu1, nu2=nu2, nu3=nu3, m12=0, m21=0, m23=m1, m32=m2, m13=0, m31=0)

    fs = Spectrum.from_phi(phi, ns, (xx,xx,xx))
    return fs

def model30_ancmig_6(params, ns, pts): # OK
    """
    Model with split between pop 1 and (2,3), with gene flow. which
    then ceases prior to 2-3 split. Split between pops 2 and 3 with gene flow, 
    ancient for-eco divergence, cessation pre forest split, but between forest block gene flow
    continues until it ceases at time 3
    for-eco symmetric, for-for asym
    """
    #10 parameters
    nu1, nuA, nu2, nu3, mA, m1, m2, T1a, T1b, T2, T3  = params

    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)


    phi = Integration.two_pops(phi, xx, T1a, nu1=nu1, nu2=nuA, m12=mA, m21=mA)

    phi = Integration.two_pops(phi, xx, T1b, nu1=nu1, nu2=nuA, m12=0, m21=0)

    phi = PhiManip.phi_2D_to_3D_split_2(xx, phi)

    phi = Integration.three_pops(phi, xx, T2, nu1=nu1, nu2=nu2, nu3=nu3, m12=0, m21=0, m23=m1, m32=m2, m13=0, m31=0)
    phi = Integration.three_pops(phi, xx, T3, nu1=nu1, nu2=nu2, nu3=nu3, m12=0, m21=0, m23=0, m32=0, m13=0, m31=0)
    fs = Spectrum.from_phi(phi, ns, (xx,xx,xx))
    return fs

def model31_ancmig_7(params, ns, pts): # OK
    """
    Model with split between pop 1 and (2,3), with gene flow. Split
    between pops 2 and 3 with gene flow, with ongoing gene flow
    between 1 and 2/3, then cessation of forest-ecotone gene flow
    with forest-forest gene flow continuing until it ceases at T4
    ancient for-eco divergence, cessation, but ancient between forest block gene flow
    for-eco sym, for-for asym
    """
    #11 parameters
    nu1, nuA, nu2, nu3, mA, m1, m2, m3, T1, T2, T3, T4 = params

    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T1, nu1=nu1, nu2=nuA, m12=mA, m21=mA)
    
    phi = PhiManip.phi_2D_to_3D_split_2(xx, phi)
    
    phi = Integration.three_pops(phi, xx, T2, nu1=nu1, nu2=nu2, nu3=nu3, m12=m1, m21=m1, m23=m2, m32=m3, m13=m3, m31=m3)

    phi = Integration.three_pops(phi, xx, T3, nu1=nu1, nu2=nu2, nu3=nu3, m12=0, m21=0, m23=m2, m32=m3, m13=0, m31=0)
    phi = Integration.three_pops(phi, xx, T4, nu1=nu1, nu2=nu2, nu3=nu3, m12=0, m21=0, m23=0, m32=0, m13=0, m31=0)
    fs = Spectrum.from_phi(phi, ns, (xx,xx,xx))
    return fs

## sym mig with size change
def model32_split_symmig_all_size(params, ns, pts): # OK
    """
    Model with split between pop 1 and (2,3), then split between 2 and 3.
    Migration is symmetrical between all population pairs (ie 1<->2, 2<->3, and 1<->3).
    with size change.
    """
    #13 parameters
    nu1, nuA, nu2, nu3, nu1b,nu2b,nu3b, mA, m1, m2, m3, T1, T2, T3 = params

    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T1, nu1=nu1, nu2=nuA, m12=mA, m21=mA)

    phi = PhiManip.phi_2D_to_3D_split_2(xx, phi)
    
    phi = Integration.three_pops(phi, xx, T2, nu1=nu1, nu2=nu2, nu3=nu3, m12=m1, m21=m1, m23=m2, m32=m2, m13=m3, m31=m3)
    phi = Integration.three_pops(phi, xx, T3, nu1=nu1b, nu2=nu2b, nu3=nu3b, m12=m1, m21=m1, m23=m2, m32=m2, m13=m3, m31=m3)
    fs = Spectrum.from_phi(phi, ns, (xx,xx,xx))
    return fs

## asym mig with size change
def model33_split_assymmig_all_size(params, ns, pts): # OK
    """
    Model with split between pop 1 and (2,3), then split between 2 and 3.
    Migration is symmetrical between all population pairs (ie 1<->2, 2<->3, and 1<->3).
    with size change.
    """
    #13 parameters
    nu1, nuA, nu2, nu3, nu1b,nu2b,nu3b, mA, m1, m2, m3, m4, T1, T2, T3 = params

    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T1, nu1=nu1, nu2=nuA, m12=mA, m21=mA)

    phi = PhiManip.phi_2D_to_3D_split_2(xx, phi)

    phi = Integration.three_pops(phi, xx, T2, nu1=nu1, nu2=nu2, nu3=nu3, m12=m1, m21=m1, m23=m2, m32=m3, m13=m4, m31=m4)
    phi = Integration.three_pops(phi, xx, T3, nu1=nu1b, nu2=nu2b, nu3=nu3b, m12=m1, m21=m1, m23=m2, m32=m3, m13=m4, m31=m4)
    fs = Spectrum.from_phi(phi, ns, (xx,xx,xx))
    return fs

#anc mig with size change
def model34_ancmig_3_size(params, ns, pts): 
    """
    Model with split between pop 1 and (2,3), with gene flow, which then stops. Split 
    between pops 2 and 3, gene flow does not occur at all.
    longest isolation, with size change.
    """
    #8 parameters
    nu1, nuA, nu2, nu3, nu1b, nu2b, nu3b, mA, T1a, T1b, T2, T3 = params

    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T1a, nu1=nu1, nu2=nuA, m12=mA, m21=mA)

    phi = Integration.two_pops(phi, xx, T1b, nu1=nu1, nu2=nuA, m12=0, m21=0)

    phi = PhiManip.phi_2D_to_3D_split_2(xx, phi)

    phi = Integration.three_pops(phi, xx, T2, nu1=nu1, nu2=nu2, nu3=nu3, m12=0, m21=0, m23=0, m32=0, m13=0, m31=0)
    phi = Integration.three_pops(phi, xx, T3, nu1=nu1b, nu2=nu2b, nu3=nu3b, m12=0, m21=0, m23=0, m32=0, m13=0, m31=0)
    fs = Spectrum.from_phi(phi, ns, (xx,xx,xx))
    return fs

def model35_ancmig_2_size(params, ns, pts): # OK
    """
    Model with split between pop 1 and (2,3), with gene flow. Split 
    between pops 2 and 3, and all gene flow ceases.
    shorter isolation
    """
    #7 parameters
    nu1, nuA, nu2, nu3, nu1b, nu2b, nu3b, mA, T1, T2, T3 = params

    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T1, nu1=nu1, nu2=nuA, m12=mA, m21=mA)

    phi = PhiManip.phi_2D_to_3D_split_2(xx, phi)
    
    phi = Integration.three_pops(phi, xx, T2, nu1=nu1, nu2=nu2, nu3=nu3, m12=0, m21=0, m23=0, m32=0, m13=0, m31=0)
    phi = Integration.three_pops(phi, xx, T3, nu1=nu1b, nu2=nu2b, nu3=nu3b, m12=0, m21=0, m23=0, m32=0, m13=0, m31=0)
    fs = Spectrum.from_phi(phi, ns, (xx,xx,xx))
    return fs

def model36_ancmig_1_size(params, ns, pts): # OK
    """
    Model with split between pop 1 and (2,3), with gene flow. Split
    between pops 2 and 3 with gene flow, then all gene flow ceases.
    shortest isolation, MODIFIED SO GENE FLOW BETEWEEN 1 and 3
    for-eco sym, for-for asym
    """
    #11 parameters
    nu1, nuA, nu2, nu3, nu1b, nu2b, nu3b, mA, m1, m2, m3, m4, T1, T2, T3, T4 = params

    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T1, nu1=nu1, nu2=nuA, m12=mA, m21=mA)
    
    phi = PhiManip.phi_2D_to_3D_split_2(xx, phi)
    
    phi = Integration.three_pops(phi, xx, T2, nu1=nu1, nu2=nu2, nu3=nu3, m12=m1, m21=m1, m23=m2, m32=m3, m13=m4, m31=m4)

    phi = Integration.three_pops(phi, xx, T3, nu1=nu1, nu2=nu2, nu3=nu3, m12=0, m21=0, m23=0, m32=0, m13=0, m31=0)
    phi = Integration.three_pops(phi, xx, T4, nu1=nu1b, nu2=nu2b, nu3=nu3b, m12=0, m21=0, m23=0, m32=0, m13=0, m31=0)
    fs = Spectrum.from_phi(phi, ns, (xx,xx,xx))
    return fs

def model37_ancmig_4_size(params, ns, pts): # OK
    """
    Model with split between pop 1 and (2,3), with gene flow. Split
    between pops 2 and 3 with gene flow, with ongoing gene flow
    between 1 and 2/3, then cessation of forest-ecotone gene flow
    ancient for-eco divergence, cessation, but between forest block gene flow
    for-eco sym, for-for asym
    """
    #16 parameters
    nu1, nuA, nu2, nu3, nu1b, nu2b, nu3b, mA, m1, m2, m3, m4, T1, T2, T3, T4 = params

    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T1, nu1=nu1, nu2=nuA, m12=mA, m21=mA)
    
    phi = PhiManip.phi_2D_to_3D_split_2(xx, phi)
    
    phi = Integration.three_pops(phi, xx, T2, nu1=nu1, nu2=nu2, nu3=nu3, m12=m1, m21=m1, m23=m2, m32=m3, m13=m4, m31=m4)

    phi = Integration.three_pops(phi, xx, T3, nu1=nu1, nu2=nu2, nu3=nu3, m12=0, m21=0, m23=m2, m32=m2, m13=0, m31=0)
    phi = Integration.three_pops(phi, xx, T4, nu1=nu1b, nu2=nu2b, nu3=nu3b, m12=0, m21=0, m23=m2, m32=m2, m13=0, m31=0)
    fs = Spectrum.from_phi(phi, ns, (xx,xx,xx))
    return fs

def model38_ancmig_5_size(params, ns, pts): # OK
    """
    Model with split between pop 1 and (2,3), with gene flow. which
    then ceases prior to 2-3 split. Split between pops 2 and 3 with gene flow, 
    ancient for-eco divergence, cessation pre forest split, but between forest block gene flow
    for-eco sym, for-for asym
    """
    #9 parameters
    nu1, nuA, nu2, nu3, nu1b, nu2b, nu3b, mA, m1, m2, T1a, T1b, T2, T3  = params

    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)


    phi = Integration.two_pops(phi, xx, T1a, nu1=nu1, nu2=nuA, m12=mA, m21=mA)

    phi = Integration.two_pops(phi, xx, T1b, nu1=nu1, nu2=nuA, m12=0, m21=0)

    phi = PhiManip.phi_2D_to_3D_split_2(xx, phi)

    phi = Integration.three_pops(phi, xx, T2, nu1=nu1, nu2=nu2, nu3=nu3, m12=0, m21=0, m23=m1, m32=m2, m13=0, m31=0)
    phi = Integration.three_pops(phi, xx, T3, nu1=nu1b, nu2=nu2b, nu3=nu3b, m12=0, m21=0, m23=m1, m32=m2, m13=0, m31=0)
    fs = Spctrum.from_phi(phi, ns, (xx,xx,xx))
    return fs


def model39_ancmig_6_size(params, ns, pts): # OK
    """
    Model with split between pop 1 and (2,3), with gene flow. which
    then ceases prior to 2-3 split. Split between pops 2 and 3 with gene flow, 
    ancient for-eco divergence, cessation pre forest split, but between forest block gene flow
    continues until it ceases at time 3
    eco-for sym, for-for asym
    """
    #10 parameters
    nu1, nuA, nu2, nu3, nu1b, nu2b, nu3b, mA, m1, m2, T1a, T1b, T2, T3, T4  = params

    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)


    phi = Integration.two_pops(phi, xx, T1a, nu1=nu1, nu2=nuA, m12=mA, m21=mA)

    phi = Integration.two_pops(phi, xx, T1b, nu1=nu1, nu2=nuA, m12=0, m21=0)

    phi = PhiManip.phi_2D_to_3D_split_2(xx, phi)

    phi = Integration.three_pops(phi, xx, T2, nu1=nu1, nu2=nu2, nu3=nu3, m12=0, m21=0, m23=m1, m32=m2, m13=0, m31=0)
    phi = Integration.three_pops(phi, xx, T3, nu1=nu1, nu2=nu2, nu3=nu3, m12=0, m21=0, m23=0, m32=0, m13=0, m31=0)
    phi = Integration.three_pops(phi, xx, T4, nu1=nu1b, nu2=nu2b, nu3=nu3b, m12=0, m21=0, m23=0, m32=0, m13=0, m31=0)
    fs = Spectrum.from_phi(phi, ns, (xx,xx,xx))
    return fs

def model40_ancmig_7_size(params, ns, pts): # OK
    """
    Model with split between pop 1 and (2,3), with gene flow. Split
    between pops 2 and 3 with gene flow, with ongoing gene flow
    between 1 and 2/3, then cessation of forest-ecotone gene flow
    with forest-forest gene flow continuing until it ceases at T4
    ancient for-eco divergence, cessation, but ancient between forest block gene flow
    for-eco sym, for-for asym
    """
    #11 parameters
    nu1, nuA, nu2, nu3, nu1b, nu2b, nu3b, mA, m1, m2, m3, m4, T1, T2, T3, T4, T5 = params

    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T1, nu1=nu1, nu2=nuA, m12=mA, m21=mA)

    phi = PhiManip.phi_2D_to_3D_split_2(xx, phi)

    phi = Integration.three_pops(phi, xx, T2, nu1=nu1, nu2=nu2, nu3=nu3, m12=m1, m21=m1, m23=m2, m32=m3, m13=4, m31=m4)

    phi = Integration.three_pops(phi, xx, T3, nu1=nu1, nu2=nu2, nu3=nu3, m12=0, m21=0, m23=m2, m32=m3, m13=0, m31=0)
    phi = Integration.three_pops(phi, xx, T4, nu1=nu1, nu2=nu2, nu3=nu3, m12=0, m21=0, m23=0, m32=0, m13=0, m31=0)
    phi = Integration.three_pops(phi, xx, T5, nu1=nu1b, nu2=nu2b, nu3=nu3b, m12=0, m21=0, m23=0, m32=0, m13=0, m31=0)
    fs = Spectrum.from_phi(phi, ns, (xx,xx,xx))
    return fs

