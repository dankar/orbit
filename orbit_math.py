import math

def degree_to_rad(deg):
    return deg * math.pi / 180
    
def rad_to_degree(rad):
    return rad * 180 / math.pi
    
def ydhms_to_s(year, day, hour, minute, seconds):
    hours = ((year-1) * 426.08 + (day-1)) * 6 + hour
    
    minute += hours * 60
    
    seconds += minute * 60
    return seconds