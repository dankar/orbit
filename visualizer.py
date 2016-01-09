import orbit
from orbit_math import *
import sys
import os
import ctypes
os.environ["PYSDL2_DLL_PATH"] = "e:\\code\\orbit\\"
from sdl2 import *

kerbol_sgp = 1172332800000000000

moho = orbit.orbit(kerbol_sgp, "Moho", 2.5263617e21, 250000.0, 5263138304.0, 0.2, 7.0, 15.0, 70.0, 3.14)
eve = orbit.orbit(kerbol_sgp, "Eve", 1.2244127e23, 700000.0, 9832684544.0, 0.01, 2.1, 0.0, 15.0, 3.14)
kerbin = orbit.orbit(kerbol_sgp, "Kerbin", 5.2915793e22, 600000.0, 13599840256.0, 0.0, 0.0, 0.0, 0.0, 3.14)
duna = orbit.orbit(kerbol_sgp, "Duna", 4.5154812e21, 320000.0, 20726155264.0, 0.051, 0.06, 0.0, 135.5, 3.14)
dres = orbit.orbit(kerbol_sgp, "Dres", 3.2191322e20, 138000.0, 40839348203.0, 0.145, 5.0, 90.0, 280.0, 3.14)
jool = orbit.orbit(kerbol_sgp, "Jool", 4.2332635e24, 6000000.0, 68773560320.0, 0.05, 1.304, 0.0, 52.0, 0.1)
eeloo = orbit.orbit(kerbol_sgp, "Eeloo", 1.1149358e21, 210000.0, 90118820000.0, 0.26, 6.15, 260.0, 50.0, 3.14)

def blit_pixel(surface, x, y, r, g, b):
    surface[(y)*(1024*4) + (x*4)] = chr(b)
    surface[(y)*(1024*4) + (x*4)+1] = chr(g)
    surface[(y)*(1024*4) + (x*4)+2] = chr(r)

def main():
    SDL_Init(SDL_INIT_VIDEO)
    window = SDL_CreateWindow(b"Orbit Visualizer",
                              SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
                              1024, 768, SDL_WINDOW_SHOWN)
    windowsurface = SDL_GetWindowSurface(window)
    
    renderer = SDL_CreateRenderer(window, -1, 0)
    
    texture = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_ARGB8888, SDL_TEXTUREACCESS_STREAMING, 1024,768);
                               
    SDL_UpdateWindowSurface(window)
    
    departure = orbit.ydhms_to_s(1, 83, 0, 0, 0)
    arrival = orbit.ydhms_to_s(1, 264, 0, 0, 0)

    g = orbit.gauss_solver(kerbin, duna, departure, arrival)
    new_orbit = g.solve()
    print "Transfer angle:", rad_to_degree(g.dv), "degrees"
    print "Delta V:", g.delta_v, "m/s"
    
    screen_buffer = ctypes.create_string_buffer('\x00' * 1024*768*4)
    
    blit_pixel(screen_buffer, 1024/2 + 0, 768/2 + 0, 255, 255, 255)
    blit_pixel(screen_buffer, 1024/2 + 1, 768/2 + 0, 255, 255, 255)
    blit_pixel(screen_buffer, 1024/2 - 1, 768/2 + 0, 255, 255, 255)
    blit_pixel(screen_buffer, 1024/2 + 0, 768/2 - 1, 255, 255, 255)
    blit_pixel(screen_buffer, 1024/2 + 0, 768/2 + 1, 255, 255, 255)
    
    print "Len:", len(screen_buffer)

    running = True
    event = SDL_Event()
    t = 0
    r = False
    while running:
        #ctypes.memset(screen_buffer, 0, 1024*768*4)
        
        v = kerbin.get_absolute_position(departure+t) / 80000000.0;
        
        blit_pixel(screen_buffer, 1024/2 + int(v.x), 768/2 + int(v.y), 0, 255, 0)
        
        v = duna.get_absolute_position(departure+t) / 80000000.0;
        
        blit_pixel(screen_buffer, 1024/2 + int(v.x), 768/2 + int(v.y), 255, 0, 0)
        
        v = new_orbit.get_absolute_position(t) / 80000000.0;
        
        blit_pixel(screen_buffer, 1024/2 + int(v.x), 768/2 + int(v.y), 255, 255, 255)
        #print v.x, v.y
        SDL_UpdateTexture(texture, None, screen_buffer, 1024 * 4);
        SDL_RenderClear(renderer);
        SDL_RenderCopy(renderer, texture, None, None)
        SDL_RenderPresent(renderer);
        
        if r and (departure + t) < arrival:
            t = t + 300
        while SDL_PollEvent(ctypes.byref(event)) != 0:
            if event.type == SDL_QUIT:
                running = False
                break
            if event.type == SDL_KEYDOWN:
                r = True
                

    SDL_DestroyWindow(window)
    SDL_Quit()
    return 0

if __name__ == "__main__":
    sys.exit(main())



#(time, deltav, angle) = find_hohmann_transfer(kerbin, duna)
#print "Hohmann solution:"
#print "Transfer time:", time, "seconds"
#print "Excess hyperbolic velocity:", deltav, "m/s"
#print "Phase angle:", rad_to_degree(angle), "degrees"

