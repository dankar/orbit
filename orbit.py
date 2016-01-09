import math
from orbit_math import *
from vector_math import *

class orbit:
    def __init__(self, sgp, _name, _mass, _radius, _semi_major_axis, _eccentricity, _inclination, _argument_of_periapsis, _longitude_of_an, _mean_anomaly):
        self.sgp = sgp
        self.name = _name
        self.mass = _mass
        self.radius = _radius
        self.semi_major_axis = _semi_major_axis
        self.eccentricity = _eccentricity
        self.inclination = degree_to_rad(_inclination)
        self.argument_of_periapsis = degree_to_rad(_argument_of_periapsis)
        self.mean_anomaly = _mean_anomaly
        self.longitude_of_an = degree_to_rad(_longitude_of_an)
        self.period = self.calc_period()
        
    def calc_period(self):
        return 2*math.pi*math.sqrt(pow(self.semi_major_axis, 3)/self.sgp)
        
    def get_future_mean_anomaly(self, seconds_from_epoch):
        n = math.sqrt(self.sgp/math.pow(self.semi_major_axis, 3))
        return self.mean_anomaly + n*seconds_from_epoch
        
    def add_seconds(self, seconds_from_epoch):
        return self.mean_anomaly + seconds_from_epoch * math.pi * 2 / self.period
        
    def get_eccentricity_anomaly(self, mean_anomaly):
        max_iterations = 30
        i = 0
        mean_anomaly = mean_anomaly / (2*math.pi)
        mean_anomaly = 2*math.pi*(mean_anomaly - math.floor(mean_anomaly))
        
        if self.eccentricity < 0.99:
            E = mean_anomaly
        else:
            E = math.pi
            
        F = E - self.eccentricity * math.sin(mean_anomaly) - mean_anomaly
        
        while math.fabs(F) > 0.00000001 and i < max_iterations:
            E = E - F / (1.0-self.eccentricity*math.cos(E))
            F = E - self.eccentricity * math.sin(E) - mean_anomaly
            i = i + 1
        return E
     
    def get_true_anomaly(self, eccentricity_anomaly):
        S = math.sin(eccentricity_anomaly)
        C = math.cos(eccentricity_anomaly)
        
        fak = math.sqrt(1.0-self.eccentricity*self.eccentricity)
        
        return math.atan2(fak*S,C-self.eccentricity)
        
     
    def get_position_vector(self, true_anomaly):
        return (self.semi_major_axis * (1-math.pow(self.eccentricity,2))) / (1+self.eccentricity*math.cos(true_anomaly))
     
    def get_flightpath_angle(self, true_anomaly):
        return math.atan2((1+self.eccentricity*math.cos(true_anomaly)), (self.eccentricity*math.sin(true_anomaly)))
        
    def get_speed(self, position_vector):
        return math.sqrt(self.sgp*(2/position_vector-1/self.semi_major_axis))
        
    def get_angular_velocity(self):
        return (2*math.pi)/self.period
        
    def transform_position(self, vec):
        lan = -self.longitude_of_an
        aop = -self.argument_of_periapsis
        inc = -self.inclination
        
        inc_axis = vector(1.0, 0.0, 0.0)
        
        vec = rotate_vector(vec, 0.0, 0.0, 1.0, lan)
        inc_axis = rotate_vector(inc_axis, 0.0, 0.0, 1.0, lan)
        vec = rotate_vector(vec, 0.0, 0.0, 1.0, aop)
        vec = rotate_vector(vec, inc_axis.x, inc_axis.y, inc_axis.z, inc)
        
        return vec
    def transform_true_anomaly(self, true_anomaly):
        return true_anomaly - self.longitude_of_an - self.argument_of_periapsis
        
    def get_absolute_position(self, seconds_from_epoch):
        mean_anomaly = -self.add_seconds(seconds_from_epoch)
        E = self.get_eccentricity_anomaly(mean_anomaly)
        true_anomaly = self.get_true_anomaly(E)
        pos = self.get_position_vector(true_anomaly)
        
        v = vector_from_angle(pos, true_anomaly)
        v = self.transform_position(v)
        return v 
        
def orbit_from_vectors(sgp, r_v, v_v, t0):
    r = r_v.magnitude()
    v = v_v.magnitude()
    
    h_v = r_v.cross(v_v)
    
    n_v = vector(-h_v.y, h_v.x, 0.0)
    
    h = h_v.magnitude()
    n = n_v.magnitude()

    e_v = (r_v * (math.pow(v, 2) - sgp / r) - (v_v * r_v.dot(v_v))) / sgp
    e = e_v.magnitude()
 
    semi_major_axis = 1.0 / (2.0 / r - math.pow(v, 2) / sgp)
    
    if semi_major_axis < 0.0:
        print "Transfer orbit is hyperbolic, orbit is currently unsupported"
        return None
    
    inclination = math.acos(h_v.z / h)
    longitude_of_an = math.acos(n_v.x / n)
    
    if n_v.y < 0:
        longitude_of_an = 2.0*math.pi - longitude_of_an

    argument_of_periapsis = 2.0*math.pi - math.acos(n_v.dot(e_v) / (n*e))

    if e_v.z < 0:
        argument_of_periapsis = 2.0*math.pi - argument_of_periapsis

    true_anomaly = math.acos((e_v.dot(r_v)) / (e*r))
    
    if r_v.dot(v_v) < 0:
        true_anomaly = 2.0*math.pi - true_anomaly
        
    cos_E = (e + math.cos(true_anomaly))/(1.0+e*math.cos(true_anomaly))
    sin_E = (math.sqrt(1.0 - math.pow(e, 2)) * math.sin(true_anomaly)) / (1.0 + e * math.cos(true_anomaly))
    
    E = math.atan2(sin_E, cos_E)
    
    mean_anomaly = E - e * math.sin(E)
    
    return orbit(sgp, "test", 0, 0, semi_major_axis, e, 180-rad_to_degree(inclination), 
                -rad_to_degree(argument_of_periapsis), -rad_to_degree(longitude_of_an), mean_anomaly)
        
class gauss_solver:
    def __init__(self, _orbit1, _orbit2, _departure_time, _arrival_time):
        self.orbit1 = _orbit1
        self.orbit2 = _orbit2
        self.departure_time = _departure_time
        self.arrival_time = _arrival_time
        self.solved = False
        
        self.true_anomaly1 = self.orbit1.get_true_anomaly(self.orbit1.get_eccentricity_anomaly(-self.orbit1.add_seconds(self.departure_time)))
        self.true_anomaly2 = self.orbit2.get_true_anomaly(self.orbit2.get_eccentricity_anomaly(-self.orbit2.add_seconds(self.arrival_time)))

        self.r1 = self.orbit1.get_position_vector(self.true_anomaly1)
        self.r2 = self.orbit2.get_position_vector(self.true_anomaly2)

        self.r1_vec = vector_from_angle(self.r1, self.true_anomaly1)
        self.r2_vec = vector_from_angle(self.r2, self.true_anomaly2)

        self.r1_vec = self.orbit1.transform_position(self.r1_vec);
        self.r2_vec = self.orbit2.transform_position(self.r2_vec);
        
        self.transform_true_anomaly1 = self.orbit1.transform_true_anomaly(self.true_anomaly1)
        self.transform_true_anomaly2 = self.orbit2.transform_true_anomaly(self.true_anomaly2)
        
        self.dv = self.transform_true_anomaly1 - self.transform_true_anomaly2
        
        #print "Transfer angle:", rad_to_degree(dv)
        
        self.k = self.calc_k()
        self.l = self.calc_l()
        self.m = self.calc_m()
        self.pi = self.calc_pi()
        self.pii = self.calc_pii()

    def calc_k(self):
        return self.r1 * self.r2 * (1-math.cos(self.dv))
        
    def calc_l(self):
        return self.r1 + self.r2
        
    def calc_m(self):
        return self.r1 * self.r2 * (1+math.cos(self.dv))
        
    def calc_pi(self):
        return self.k / (self.l + math.sqrt(2*self.m))
        
    def calc_pii(self):
        return self.k / (self.l - math.sqrt(2*self.m))
        
    def calc_a(self, p):
        return self.m*self.k*p/((2.0*self.m - math.pow(self.l, 2)) * math.pow(p, 2) + 2.0 * self.k * self.l * p - math.pow(self.k, 2))
        
    def calc_f(self, p):
        return 1.0 - self.r2 / p * (1.0 - math.cos(self.dv))
        
    def calc_fdot(self, p):
        return math.sqrt(self.orbit1.sgp/p) * math.tan(self.dv/2) * ((1-math.cos(self.dv))/p - 1.0/self.r1 - 1/self.r2)
        
    def calc_g(self, p):
        return self.r1 * self.r2 * math.sin(self.dv) / math.sqrt(self.orbit1.sgp * p)
        
    def calc_dF(self, a, f):
        return math.acosh(1.0 - self.r1 / a * (1.0 - f))
        
    def calc_cos_dE(self, a, f):
        return 1.0 - (self.r1 / a * (1.0 - f))
        
    def calc_sin_dE(self, a, fdot):
        return -self.r1 * self.r2 * fdot / math.sqrt(self.orbit1.sgp*a)
        
    def calc_t_from_dE(self, a, dE, g):
        return g + math.sqrt(math.pow(a, 3) / self.orbit1.sgp) * (dE - math.sin(dE))
        
    def calc_t_from_dF(self, a, dF, g):
        return g + math.sqrt(math.pow(-a, 3)/self.orbit1.sgp) * (math.sinh(dF) - dF)
        
    def calc_at(self, p):
        a = self.calc_a(p)
        f = self.calc_f(p)
        fdot = self.calc_fdot(p)
        g = self.calc_g(p)

        if a < 0.0:
            dF = self.calc_dF(a, f)
            t = self.calc_t_from_dF(a, dF, g)
        else:
            cosdE = self.calc_cos_dE(a, f)
            sindE = self.calc_sin_dE(a, fdot)
            
            dE = math.atan2(sindE, cosdE)
            
            t = self.calc_t_from_dE(a, dE, g)

        return (a, t)
        
    def calc_new_p(self, p1, p0, t_error1, t_error0):
        delta_p = -((p0 - p1) * t_error1)/(t_error0 - t_error1)
        
        step_limiter = math.pow(delta_p / p1, 2) + 1.0
        
        return p1 + delta_p / step_limiter
        
    def solve(self):
        p0 = (self.pi*2.0+self.pii)/3.0
        p1 = (self.pii*2.0+self.pi)/3.0
        
        max_iterations = 50
        iterations = 0
        
        (a, t0) = self.calc_at(p0)
        
        while True:
            (a, t1) = self.calc_at(p1)
            
            t_error0 = t0 - (self.arrival_time - self.departure_time)
            t_error1 = t1 - (self.arrival_time - self.departure_time)
            
            p = self.calc_new_p(p1, p0, t_error1, t_error0) 
            
            if math.fabs(t_error1) < 1.0:
                break
                
            p0 = p1
            p1 = p
            t0 = t1
            
            iterations = iterations + 1
            
            if iterations > max_iterations:
                print "No solution"
                self.solved = 0
                self.delta_v = 0.0
                return 0
        
        f = self.calc_f(p1)
        g = self.calc_g(p1)
        
        v1 = (self.r2_vec - (self.r1_vec * f)) / g
        
        fpa = self.transform_true_anomaly1 - self.orbit1.get_flightpath_angle(self.true_anomaly1)
        pos = self.orbit1.get_position_vector(self.true_anomaly1)
        speed = self.orbit1.get_speed(pos)
        planet_v = vector_from_angle(speed, fpa)
        
        o = orbit_from_vectors(self.orbit1.sgp, self.r1_vec, v1, self.departure_time)
        #o = orbit_from_vectors(1.327124e20, vector(7.079944e10, -1.345206e11, 0.0), vector(28996.2, 15232.7, 1289.2))
        v1 = v1 - planet_v
        self.delta_v = v1.magnitude()
        return o

def find_hohmann_transfer(orbit1, orbit2):
    angular_vel = orbit2.get_angular_velocity()
    travel_time = math.pi * math.sqrt(math.pow(orbit1.semi_major_axis+orbit2.semi_major_axis, 3) / (8*orbit1.sgp))
    excess_v = math.sqrt(orbit1.sgp/orbit1.semi_major_axis) * (math.sqrt((2.0*orbit2.semi_major_axis)/(orbit1.semi_major_axis+orbit2.semi_major_axis))-1)
    phase_angle = travel_time * angular_vel

    return (travel_time, excess_v, math.pi - phase_angle)
