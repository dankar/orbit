import math

class vector:
    def __init__(self, _x, _y, _z):
        self.x = _x
        self.y = _y
        self.z = _z
        
    def magnitude(self):
        return math.sqrt(self.x*self.x + self.y*self.y + self.z*self.z)
        
    def __mul__(self, scalar):
        return vector(self.x * scalar, self.y * scalar, self.z * scalar)
        
    def __div__(self, scalar):
        return vector(self.x / scalar, self.y / scalar, self.z / scalar) 
        
    def __sub__(self, other):
        return vector(self.x - other.x, self.y - other.y, self.z - other.z)
        
    def cross(self, other):
        return vector(self.y * other.z - self.z * other.y, self.z * other.x - self.x * other.z, self.x * other.y - self.y * other.x)
        
    def dot(self, other):
        return self.x * other.x + self.y * other.y + self.z * other.z
        
def vector_from_angle(pos, angle):
        return vector(math.cos(angle) * pos, math.sin(angle) * pos, 0.0)
        
class quaternion:
    def __init__(self, _w, _x, _y, _z):
        self.w = _w
        self.x = _x
        self.y = _y
        self.z = _z
        
    def conjugate(self):
        return quaternion(self.w, -self.x, -self.y, -self.z)
        
    def magnitude(self):
        return math.sqrt(self.x*self.x + self.y*self.y + self.z*self.z)
        
    def __div__(self, scalar):
        return quaternion(self.w/scalar, self.x/scalar, self.y/scalar, self.z/scalar)
        
    def normalized(self):
        mag = self.magnitude();
        return self/mag;
        
    def __mul__(self, other):
        return quaternion(self.w * other.w - self.x * other.x - self.y * other.y - self.z * other.z,
                        self.w * other.x + self.x * other.w + self.y * other.z - self.z * other.y,
                        self.w * other.y + self.y * other.w + self.z * other.x - self.x * other.z,
                        self.w * other.z + self.z * other.w + self.x * other.y - self.y * other.x)

def quaternion_from_rotation(x, y, z, angle):
    rot = quaternion(0.0, x, y, z).normalized()
    angle = angle / 2
    rot.w = math.cos(angle)
    rot.x = rot.x * math.sin(angle)
    rot.y = rot.y * math.sin(angle)
    rot.z = rot.z * math.sin(angle)
    
    return rot

def rotate_vector(vec, x, y, z, angle):
    q = quaternion_from_rotation(x, y, z, angle)
    
    q_conj = q.conjugate()
    qvec = quaternion(0.0, vec.x, vec.y, vec.z)
    result = q * qvec * q_conj

    return vector(result.x, result.y, result.z)
    
    