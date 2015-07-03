import sys
import numpy
from galpy.potential import MWPotential2014
from galpy.orbit import Orbit
from galpy.util import bovy_conversion
def get_initial_condition(to):
    """Find an initial condition near apocenter that ends up today near the end of the stream and is close to a gyrfalcon stepsize
    to= Find an apocenter that is just more recent than this time"""
    vo, ro= 220., 8.
    # Center of the stream
    o= Orbit([0.10035165,-0.81302488,0.80986668,0.58024425,0.92753945,
               0.88763126],ro=ro,vo=vo)
    #Flip for backwards integration
    of= o.flip()
    ts= numpy.linspace(0.,to/bovy_conversion.time_in_Gyr(vo,ro),10001)
    of.integrate(ts,MWPotential2014)
    # Find the closest apocenter to the time to
    rs= numpy.sqrt(of.R(ts)**2.+of.z(ts)**2.)
    drs= rs-numpy.roll(rs,1)
    nearApo= (drs < 0.)*(numpy.roll(drs,1) > 0.)
    tsNearApo= ts[nearApo]
    tsNearApo*= bovy_conversion.time_in_Gyr(vo,ro)
    tfgf= numpy.amax(tsNearApo)
    #Round to nearest 2.**-8.
    tfgf= round(tfgf/2.**-8.)*2.**-8
    tf= tfgf/bovy_conversion.time_in_Gyr(vo,ro)
    print 'Current: %s,%s,%s,%s,%s,%s' % (of.x()[0], of.y()[0], of.z(), -of.vx()[0], -of.vy()[0], -of.vz())
    print 'At %g Gyr: %s,%s,%s,%s,%s,%s' % (tf*bovy_conversion.time_in_Gyr(vo,ro),of.x(tf)[0], of.y(tf)[0], of.z(tf), -of.vx(tf,use_physical=False)[0]*bovy_conversion.velocity_in_kpcGyr(vo,ro), -of.vy(tf,use_physical=False)[0]*bovy_conversion.velocity_in_kpcGyr(vo,ro), -of.vz(tf,use_physical=False)*bovy_conversion.velocity_in_kpcGyr(vo,ro))
    return None

if __name__ == '__main__':
    get_initial_condition(float(sys.argv[1]))
