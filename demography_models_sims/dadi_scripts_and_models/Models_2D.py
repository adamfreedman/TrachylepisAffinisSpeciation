import numpy
from dadi import Numerics, PhiManip, Integration
from dadi.Spectrum_mod import Spectrum

'''
Models for testing two population scenarios.
'''

def model20_anc_sym_mig_size_3epoch_genflow_2levels_and_secondarycontact(params, ns, pts):
    """
    Model with split and no gene flow, followed by asymmetrical gene flow, size change.

    nu1a: Size of population 1 after split.
    nu2a: Size of population 2 after split.
    T1a: Time in the past of split (in units of 2*Na generations) where gene flow ceases
    T1b: Time post-split at which gene flow resumes upon secondary contact post expansion
    nu1b: Size of population 1 after time interval.
    nu2b: Size of population 2 after time interval.
    nu1c: Size of population 1 after time interval.
    nu2c: Size of population 2 after time interval.
    T2: The scale time between the ancient migration and present.
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    """
    nu1a, nu2a, nu1b, nu2b, nu1c, nu2c, ma, mb, T1a, T1b, T2 = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T1a, nu1a, nu2a, m12=ma, m21=ma)
    phi = Integration.two_pops(phi, xx, T1b, nu1b, nu2b, m12=0, m21=0)

    phi = Integration.two_pops(phi, xx, T2, nu1c, nu2c, m12=mb, m21=mb)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

def model19_anc_asym_mig_size_3epoch_genflow_2levels_and_secondarycontact(params, ns, pts):
    """
    Model with split and no gene flow, followed by asymmetrical gene flow, size change.

    nu1a: Size of population 1 after split.
    nu2a: Size of population 2 after split.
    T1a: Time in the past of split (in units of 2*Na generations) where gene flow ceases
    T1b: Time post-split at which gene flow resumes upon secondary contact post expansion
    nu1b: Size of population 1 after time interval.
    nu2b: Size of population 2 after time interval.
    nu1c: Size of population 1 after time interval.
    nu2c: Size of population 2 after time interval.
    T2: The scale time between the ancient migration and present.
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    """
    nu1a, nu2a, nu1b, nu2b, nu1c, nu2c, m12a, m21a, m12b, m21b,T1a, T1b, T2 = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T1a, nu1a, nu2a, m12=m12a, m21=m21a)
    phi = Integration.two_pops(phi, xx, T1b, nu1b, nu2b, m12=0, m21=0)

    phi = Integration.two_pops(phi, xx, T2, nu1c, nu2c, m12=m12b, m21=m21b)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs



def model18_anc_asym_mig_size_3epoch_genflow_and_secondarycontact(params, ns, pts):
    """
    Model with split and no gene flow, followed by asymmetrical gene flow, size change.

    nu1a: Size of population 1 after split.
    nu2a: Size of population 2 after split.
    T1a: Time in the past of split (in units of 2*Na generations) where gene flow ceases
    T1b: Time post-split at which gene flow resumes upon secondary contact post expansion
    nu1b: Size of population 1 after time interval.
    nu2b: Size of population 2 after time interval.
    nu1c: Size of population 1 after time interval.
    nu2c: Size of population 2 after time interval.
    T2: The scale time between the ancient migration and present.
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    """
    nu1a, nu2a, nu1b, nu2b, nu1c, nu2c, m12, m21, T1a, T1b, T2 = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T1a, nu1a, nu2a, m12=m12, m21=m21)
    phi = Integration.two_pops(phi, xx, T1b, nu1b, nu2b, m12=0, m21=0)

    phi = Integration.two_pops(phi, xx, T2, nu1c, nu2c, m12=m12, m21=m21)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs
def model17_anc_asym_mig_size_3epoch_gradient_gflowdescent2(params, ns, pts):
    """
    Model with split and no gene flow, followed by asymmetrical gene flow, size change.

    nu1a: Size of population 1 after split.
    nu2a: Size of population 2 after split.
    T1a: Time in the past of split (in units of 2*Na generations) where gene flow ceases
    T1b: Time post-split at which gene flow resumes upon secondary contact post expansion
    nu1b: Size of population 1 after time interval.
    nu2b: Size of population 2 after time interval.
    T2: The scale time between the ancient migration and present.
    m12a: Migration from pop 2 to pop 1 (2*Na*m12).
    m21a: Migration from pop 1 to pop 2.
    m12b
    m21b
    """
    nu1a, nu2a, nu1b, nu2b, m12a, m21a, m12b, m21b, T1a, T1b, T2 = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T1a, nu1a, nu2a, m12=m12a, m21=m21a)
    phi = Integration.two_pops(phi, xx, T1b, nu1b, nu2b, m12=m12b, m21=m21b)

    phi = Integration.two_pops(phi, xx, T2, nu1b, nu2b, m12=0, m21=0)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

def model16_anc_asym_mig_size_3epoch_gradient_gflowdescent1(params, ns, pts):
    """
    Model with split and no gene flow, followed by asymmetrical gene flow, size change.

    nu1a: Size of population 1 after split.
    nu2a: Size of population 2 after split.
    T1a: Time in the past of split (in units of 2*Na generations) where gene flow ceases
    T1b: Time post-split at which gene flow resumes upon secondary contact post expansion
    nu1b: Size of population 1 after time interval.
    nu2b: Size of population 2 after time interval.
    nu1c: Size of population 1 after time interval.
    nu2c: Size of population 2 after time interval.
    T2: The scale time between the ancient migration and present.
    m12a: Migration from pop 2 to pop 1 (2*Na*m12).
    m21a: Migration from pop 1 to pop 2.
    m12b
    m21b
    """
    nu1a, nu2a, nu1b, nu2b, nu1c, nu2c, m12a, m21a, m12b, m21b, T1a, T1b, T2 = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T1a, nu1a, nu2a, m12=m12a, m21=m21a)
    phi = Integration.two_pops(phi, xx, T1b, nu1b, nu2b, m12=m12b, m21=m21b)

    phi = Integration.two_pops(phi, xx, T2, nu1c, nu2c, m12=0, m21=0)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

def model21_anc_asym_mig_size_3epoch_earlysecondarycontact(params, ns, pts):
    """
    Model with split and no gene flow, followed by asymmetrical gene flow, size change.

    nu1a: Size of population 1 after split.
    nu2a: Size of population 2 after split.
    T1a: Time in the past of split (in units of 2*Na generations) where gene flow ceases
    T1b: Time post-split at which gene flow resumes upon secondary contact post expansion
    nu1b: Size of population 1 after time interval.
    nu2b: Size of population 2 after time interval.
    nu1c: Size of population 1 after time interval.
    nu2c: Size of population 2 after time interval.
    T2: The scale time between the ancient migration and present.
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    """
    nu1a, nu2a, nu1b, nu2b, nu1c, nu2c, m12, m21, T1a, T1b, T2 = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T1a, nu1a, nu2a, m12=0, m21=0)
    phi = Integration.two_pops(phi, xx, T1b, nu1b, nu2b, m12=m12, m21=m21)

    phi = Integration.two_pops(phi, xx, T2, nu1c, nu2c, m12=0, m21=0)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs
    

def model15_anc_asym_mig_size_3epoch_gradient(params, ns, pts):
    """
    Model with split and no gene flow, followed by asymmetrical gene flow, size change.

    nu1a: Size of population 1 after split.
    nu2a: Size of population 2 after split.
    T1a: Time in the past of split (in units of 2*Na generations) where gene flow ceases
    T1b: Time post-split at which gene flow resumes upon secondary contact post expansion
    nu1b: Size of population 1 after time interval.
    nu2b: Size of population 2 after time interval.
    nu1c: Size of population 1 after time interval.
    nu2c: Size of population 2 after time interval.
    T2: The scale time between the ancient migration and present.
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    """
    nu1a, nu2a, nu1b, nu2b, nu1c, nu2c, m12, m21, T1a, T1b, T2 = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T1a, nu1a, nu2a, m12=m12, m21=m21)
    phi = Integration.two_pops(phi, xx, T1b, nu1b, nu2b, m12=0, m21=0)

    phi = Integration.two_pops(phi, xx, T2, nu1c, nu2c, m12=0, m21=0)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs
        

def no_divergence(notused, ns, pts):
    """
    Standard neutral model, populations never diverge.
    """
    
    xx = Numerics.default_grid(pts)
    
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs


def model1_no_mig(params, ns, pts):
    """
    Split into two populations, no migration.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    T: Time in the past of split (in units of 2*Na generations) 
    """
    nu1, nu2, T = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T, nu1, nu2, m12=0, m21=0)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs


def model2_sym_mig(params, ns, pts):
    """
    Split into two populations, with symmetric migration.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    T: Time in the past of split (in units of 2*Na generations) 
    m: Migration rate between populations (2*Na*m)
    """
    nu1, nu2, m, T = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T, nu1, nu2, m12=m, m21=m)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs


def model3_asym_mig(params, ns, pts):
    """
    Split into two populations, with different migration rates.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    T: Time in the past of split (in units of 2*Na generations) 
    m12: Migration from pop 2 to pop 1 (2*Na*m12)
    m21: Migration from pop 1 to pop 2
	"""
    nu1, nu2, m12, m21, T = params
    xx = Numerics.default_grid(pts)
    
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    
    phi = Integration.two_pops(phi, xx, T, nu1, nu2, m12=m12, m21=m21)
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    
    return fs    

def model4_anc_sym_mig(params, ns, pts):
    """
    Model with split and symmetric migration followed by isolation.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m: Migration between pop 2 and pop 1.
    T1: The scaled time between the split and the ancient migration (in units of 2*Na generations).
    T2: The scaled time between the ancient migration and present.
    """
    nu1, nu2, m, T1, T2 = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    
    phi = Integration.two_pops(phi, xx, T1, nu1, nu2, m12=m, m21=m)

    phi = Integration.two_pops(phi, xx, T2, nu1, nu2, m12=0, m21=0)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs
 
def model5_anc_asym_mig(params, ns, pts):
    """
    Model with split and asymmetric migration followed by isolation.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    T1: The scaled time between the split and the ancient migration (in units of 2*Na generations).
    T2: The scaled time between the ancient migration and present.
    """
    nu1, nu2, m12, m21, T1, T2 = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    
    phi = Integration.two_pops(phi, xx, T1, nu1, nu2, m12=m12, m21=m21)
    
    phi = Integration.two_pops(phi, xx, T2, nu1, nu2, m12=0, m21=0)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

def model6_sec_contact_sym_mig(params, ns, pts):
    """
    Model with split and no gene flow, followed by period of symmetrical gene flow.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m: Migration between pop 2 and pop 1.
    T1: The scaled time between the split and the secondary contact (in units of 2*Na generations).
    T2: The scaled time between the secondary contact and present.
    """
    nu1, nu2, m, T1, T2 = params

    xx = Numerics.default_grid(pts)
    
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T1, nu1, nu2, m12=0, m21=0)

    phi = Integration.two_pops(phi, xx, T2, nu1, nu2, m12=m, m21=m)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

def model7_sec_contact_asym_mig(params, ns, pts):
    """
    Model with split and no gene flow, followed by period of asymmetrical gene flow.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    T1: The scaled time between the split and the secondary contact (in units of 2*Na generations).
    T2: The scaled time between the secondary contact and present.
    """
    nu1, nu2, m12, m21, T1, T2 = params

    xx = Numerics.default_grid(pts)
    
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T1, nu1, nu2, m12=0, m21=0)

    phi = Integration.two_pops(phi, xx, T2, nu1, nu2, m12=m12, m21=m21)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs



#######################################################################################################
#Models involving size changes

def model8_no_mig_size(params, ns, pts):
    """
    Split into two populations with no migration, then size change.

    nu1a: Size of population 1 after split.
    nu2a: Size of population 2 after split.
    T1: Time in the past of split (in units of 2*Na generations)
    nu1b: Size of population 1 after time interval.
    nu2b: Size of population 2 after time interval.
    T2: Time of population size change.
    """
    nu1a, nu2a, nu1b, nu2b, T1, T2 = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T1, nu1a, nu2a, m12=0, m21=0)
    
    phi = Integration.two_pops(phi, xx, T2, nu1b, nu2b, m12=0, m21=0)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

def model9_sym_mig_size(params, ns, pts):
    """
    Split into two populations with symmetric migration, size change.

    nu1a: Size of population 1 after split.
    nu2a: Size of population 2 after split.
    T1: Time in the past of split (in units of 2*Na generations)
    nu1b: Size of population 1 after time interval.
    nu2b: Size of population 2 after time interval.
    T2: Time of population size change.
    m: Migration rate between populations (2*Na*m)
    """
    nu1a, nu2a, nu1b, nu2b, m, T1, T2 = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T1, nu1a, nu2a, m12=m, m21=m)
    
    phi = Integration.two_pops(phi, xx, T2, nu1b, nu2b, m12=m, m21=m)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

def model10_asym_mig_size(params, ns, pts):
    """
    Split into two populations with different migration rates, size change.

    nu1a: Size of population 1 after split.
    nu2a: Size of population 2 after split.
    T1: Time in the past of split (in units of 2*Na generations)
    nu1b: Size of population 1 after time interval.
    nu2b: Size of population 2 after time interval.
    T2: Time of population size change.
    m12: Migration from pop 2 to pop 1 (2*Na*m12)
    m21: Migration from pop 1 to pop 2
	"""
    nu1a, nu2a, nu1b, nu2b, m12, m21, T1, T2 = params
    xx = Numerics.default_grid(pts)
    
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    
    phi = Integration.two_pops(phi, xx, T1, nu1a, nu2a, m12=m12, m21=m21)
    
    phi = Integration.two_pops(phi, xx, T2, nu1b, nu2b, m12=m12, m21=m21)
    
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    
    return fs    

def model11_anc_sym_mig_size(params, ns, pts):
    """
    Model with split and no gene flow, followed by symmetrical gene flow, size change.

    nu1a: Size of population 1 after split.
    nu2a: Size of population 2 after split.
    T1: Time in the past of split (in units of 2*Na generations)
    nu1b: Size of population 1 after time interval.
    nu2b: Size of population 2 after time interval.
    T2: The scale time between the ancient migration and present.
    m: Migration between pop 2 and pop 1.
    """
    nu1a, nu2a, nu1b, nu2b, m, T1, T2 = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    
    phi = Integration.two_pops(phi, xx, T1, nu1a, nu2a, m12=m, m21=m)

    phi = Integration.two_pops(phi, xx, T2, nu1b, nu2b, m12=0, m21=0)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

def model12_anc_asym_mig_size(params, ns, pts):
    """
    Model with split and no gene flow, followed by asymmetrical gene flow, size change.

    nu1a: Size of population 1 after split.
    nu2a: Size of population 2 after split.
    T1: Time in the past of split (in units of 2*Na generations)
    nu1b: Size of population 1 after time interval.
    nu2b: Size of population 2 after time interval.
    T2: The scale time between the ancient migration and present.
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    """
    nu1a, nu2a, nu1b, nu2b, m12, m21, T1, T2 = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    
    phi = Integration.two_pops(phi, xx, T1, nu1a, nu2a, m12=m12, m21=m21)
    
    phi = Integration.two_pops(phi, xx, T2, nu1b, nu2b, m12=0, m21=0)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

def model13_sec_contact_sym_mig_size(params, ns, pts):
    """
    Model with split in isolation, followed by secondary contact with asymmetrical gene flow, size change.

    nu1a: Size of population 1 after split.
    nu2a: Size of population 2 after split.
    T1: The scaled time between the split and the secondary contact (in units of 2*Na generations).
    nu1b: Size of population 1 after time interval.
    nu2b: Size of population 2 after time interval.
    T2: The scale time between the secondary contact and present.
    m: Migration between pop 2 and pop 1.
    """
    nu1a, nu2a, nu1b, nu2b, m, T1, T2 = params

    xx = Numerics.default_grid(pts)
    
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T1, nu1a, nu2a, m12=0, m21=0)

    phi = Integration.two_pops(phi, xx, T2, nu1b, nu2b, m12=m, m21=m)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs
    
def model14_sec_contact_asym_mig_size(params, ns, pts):
    """
    Model with split in isolation, followed by secondary contact with asymmetrical gene flow, size change.

    nu1a: Size of population 1 after split.
    nu2a: Size of population 2 after split.
    T1: The scaled time between the split and the secondary contact (in units of 2*Na generations).
    nu1b: Size of population 1 after time interval.
    nu2b: Size of population 2 after time interval.
    T2: The scale time between the secondary contact and present.
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    """
    nu1a, nu2a, nu1b, nu2b, m12, m21, T1, T2 = params

    xx = Numerics.default_grid(pts)
    
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    
    phi = Integration.two_pops(phi, xx, T1, nu1a, nu2a, m12=0, m21=0)

    phi = Integration.two_pops(phi, xx, T2, nu1b, nu2b, m12=m12, m21=m21)
    
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs
