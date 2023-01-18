import numpy as np
from PIL import Image, ImageDraw, ImageFont
from os import listdir, scandir, makedirs
from sys import argv

# init
i = 1
width, height = (3000, 1431) # manually define width, height
if len(argv) == 3: # overwrite width, height if cl arguments present
    width = int(argv[1])
    height = int(argv[2])
print(width, height, argv)
frame = np.empty((height*2, width*2, 3), np.uint8)
# keeps track of which videos still have more frames
unfinished = [True] * 4 
# this fontpath is for Linux, if you have font related errors changing this
# path is a good place to start
fontpath = '/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf'
font = ImageFont.truetype(fontpath, 72)
frames = np.empty((4, height, width, 3), np.uint8)

# find what values were used
subdirs = [dir.name for dir in scandir('.') if dir.is_dir() and 'frames_' in dir.name]
if len(subdirs) != 4:
    errmsg = f"The number of folders found was incorrect, please make sure there are exactly four. Found {len(subdirs)} folders at {subdirs}."
    raise ValueError(errmsg)

subscripts = [dirname.split('_')[1] for dirname in subdirs]

# sort the frames
inds = np.argsort([float(s) for s in subscripts])
subscripts = [subscripts[i] for i in inds]
subdirs = [subdirs[i] for i in inds]

numFrames = max([len(listdir(f'{dir}')) for dir in subdirs])

makedirs('frames', exist_ok=True)

while any(unfinished):
    # extract frames
    for j, dir in enumerate(subdirs):
        try: 
            with Image.open(f'{dir}/frame_fireLine_{i:04}.png') as img:
                frames[j] = np.array(img)
        except FileNotFoundError: unfinished[j] = False

    # add frames into single 4x4 frame
    frame[0:height, 0:width, :] = frames[0]
    frame[0:height, width:, :] = frames[1]
    frame[height:, 0:width, :] = frames[2]
    frame[height:, width:, :] = frames[3]

    # convert array into PIL image
    img = Image.fromarray(frame)
    # annotate 4x4 with property value
    draw = ImageDraw.Draw(img)
    draw.text((0, 0), f"slopeStrength = {subscripts[0]}", (0,0,0), font=font)
    draw.text((width, 0), f"slopeStrength = {subscripts[1]}", (0,0,0), font=font)
    draw.text((0, height), f"slopeStrength = {subscripts[2]}", (0,0,0), font=font)
    draw.text((width, height), f"slopeStrength = {subscripts[3]}", (0,0,0), font=font)

    # save image
    img.save(f'frames/fire4x4_{i:04d}.png')
    img.close()
    print(100 * i // numFrames, '%', sep='')
    i += 1
