load data/interim/refactor-test/22/H2O/1_make_topology/mc_sol.pdb

load_traj data/interim/refactor-test/22/H2O/11_GaMD_full/100/0/10236f52ef18b2a3_traj.netcdf, format=nc, start=1, stop=100

bg_color black
set ray_opaque_background, off
set orthoscopic, 0
set transparency, 0.1
set dash_gap, 0
set ray_trace_mode, 1
set ray_texture, 2
set antialias, 3
set ambient, 0.5
set spec_count, 5
set shininess, 50
set specular, 1
set reflect, .1
space cmyk



remove solvent
remove metals
hide cartoon

show sticks
show spheres
color gray85,elem C
color gray98,elem H
color slate,elem N
color red,elem O
set stick_radius,0.07
set sphere_scale,0.18
set sphere_scale,0.13, elem H

# preset.ball_and_stick(selection='all', mode=1)

orient
center

set antialias, 5 # for better ray image quality
set scene_animation_duration, 5 # the longer the animation, the smoother the movie
set ray_trace_frames, on # to render the individual images


set_view (\
    -0.854636133,   -0.333118379,   -0.398282081,\
     0.441424698,   -0.870027900,   -0.219531968,\
    -0.273384601,   -0.363431931,    0.890604794,\
     0.000000000,    0.000000000,  -37.951503754,\
     1.766034603,  -10.998698235,    9.105729103,\
    29.921253204,   45.981754303,  -20.000000000 )


mset 1x200

mview store, 1, state=1, object=mc_sol

mview store, 100, state=100, object=mc_sol

mview store, 200, state=1, object=mc_sol

intra_fit mc_sol
smooth

# mpng movie/, width=1920, height=1080
# smooth

